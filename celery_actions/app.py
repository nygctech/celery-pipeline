import sys

from flask import Flask, render_template, request, jsonify
import subprocess
import requests
import os
import re
import json
import numpy as np
import glob
import scanpy as sc
import traceback
import pandas as pd
import random
import string


app = Flask(__name__)

CLOUD_ADDRESS="https://us-central1-techinno.cloudfunctions.net/"
MAX_IMAGE_WIDTH=9500


@app.route('/')
def index():
    return render_template('index.html')

@app.route('/add_project')
def add_project():
    return render_template('add_project.html')

@app.route('/add_sample')
def add_sample():
    return render_template('add_sample.html')


@app.route('/add_sample_do', methods=['POST'])
def add_sample_do():
    #this change of cloud address was in case it ever changed, but it's unlikely, should be removed gradually
    cloud_address=CLOUD_ADDRESS
    error_messages={
        "token":  "Token is needed to authenticate request.",
        "name" :  "Project name is needed.",
        "rid"  :  "Project needs a role.",
        "classes":"Project needs Classes."
    }


def generate_random_string(length):
    # Define the characters to choose from
    characters = string.ascii_letters + string.digits
    # Generate a random string
    random_string = ''.join(random.choice(characters) for _ in range(length))
    return random_string

def parse_identify_str(ident):
    s=ident.lower()
    pattbrackets = r'\[\d+\]'
    pattsizes = r'(\d+)x(\d+)'
    level_dict={}
    if "[0]" in s:
        levels = re.finditer(pattbrackets, s)
        for match in levels:
            l=match.group()
            substr=s[match.start():]
            size=re.search(pattsizes, substr)
            if size:
                level_dict[l]=[int(g) for g in size.groups()]
                print(l, level_dict[l])

        widths=[(level_dict[l][0],l) for l in level_dict]
        windices=np.flip(np.argsort([a[0] for a in widths]))

        ld={}
        for wi in windices:
            l=widths[wi][1]
            ld[l]=level_dict[l]
        return ld
    else:
        size=re.search(pattsizes, s)
        level_dict["[unique]"]=[int(g) for g in size.groups()]
    return level_dict

def convert_cmd(json_data, levels, mainlevel, maxwidth):
    try:
        #module=json_data["module"]
        module="ImageMagick/7.1.1-15-GCCcore-12.3.0"
        imageloc=json_data["image_location"]
        loadmodcmd=f"module load {module}"
        scratchloc="/scratch"

        imagename = imageloc.split("/")[-1]
        namenoformat = "".join(imagename.split(".")[:-1]) or imagename

        is_jpeg=False
        is_width=False #is width good enough to pass
        cmd=None
        scale=1.0

        maxlevelwidth=0
        for l in levels:
            if levels[l][0]>maxlevelwidth:
                maxlevelwidth= levels[l][0]

        mainlevel_width=levels[mainlevel][0]
        if mainlevel_width<maxwidth:
            is_width=True

        entry_file_size_bytes = os.path.getsize(imageloc)
        entry_file_size_mb = entry_file_size_bytes / (1024*1024)

        if any(x in imageloc.lower() for x in [".jpg","jpeg",".jp2"]):
            is_jpeg=True

        if is_jpeg and entry_file_size_mb < 35 and is_width:
            #upload as is
            return {"stdout":"","convertloc":imageloc,"scale":scale,"message":"No need to do anything, upload as is"},None,None

        convertedname=f"{scratchloc}/{namenoformat}.jpg"
        message=""
        if mainlevel == "[unique]":
            if entry_file_size_mb < 35 and is_width:
                cmd=f"convert {imageloc} {convertedname}"
                message="just convert to jpeg"
            else:
                cmd=f"convert {imageloc} -resize {maxwidth}x {convertedname}"
                scale=maxwidth/maxlevelwidth
                message=f"convert to width {maxwidth} then to jpeg"
        else:
            #cmd=f"convert {imageloc}{mainlevel} -resize {maxwidth}x -define jpeg:extent=35500 {convertedname}"
            #the define extent brings problems sometimes, maybe just not use it
            cmd=f"convert {imageloc}{mainlevel} -resize {maxwidth}x {convertedname}"
            print("would've chosen to reduce jpeg, but forget it ")
            scale=maxwidth/maxlevelwidth
            message=f"convert level {mainlevel} to width {maxwidth} then to jpeg"



        scale=maxwidth/maxlevelwidth

        full_command = f"{loadmodcmd}"
        full_command+= f" && {cmd}"

        hostname=os.uname().nodename

        if "ne1" not in hostname:
            full_command = f"identify {imageloc}"

        responsedict={}
        process = subprocess.Popen(full_command,
                                   shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   executable='/bin/bash')

        stdout, stderr = process.communicate()

        return {"stdout":stdout.decode(),"convertloc":convertedname,"scale":scale,"message":f" {message}. Did this command: '{cmd}'"},stderr,None

    except Exception as e:
        print("Exception in convert was:",e)
        print(traceback.format_exc())
        return {"stdout":stdout},stderr,e

def identify_cmd(json_data):
    try:
        #module=json_data["module"]
        module="ImageMagick/7.1.1-15-GCCcore-12.3.0"
        imageloc=json_data["image_location"]
        loadmodcmd=f"module load {module}"

        full_command = f"{loadmodcmd} && identify {imageloc}"

        hostname=os.uname().nodename

        if "ne1" not in hostname:
            full_command = f"identify {imageloc}"

        responsedict={}
        process = subprocess.Popen(full_command,
                                   shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   executable='/bin/bash')

        stdout, stderr = process.communicate()

        return stdout,stderr,None

    except Exception as e:
        print("Exception in identify:", e)
        print(traceback.format_exc())
        return stdout,stderr,e

@app.route('/convert_image', methods=['POST'])
def convert_image():
    if not request.is_json:
        return jsonify({"error": "Request is not json formatted"}), 400

    json_data = request.json
    scale, width = None, None
    message = ""
    missing = []

    def valid_s(ese):
        if isinstance(ese,str):
            return ese.isdigit()
        elif isinstance(ese,float):
            return ese >0.0
        return False

    def valid_w(dw):
        if isinstance(dw,str):
            return dw.isdigit() and "." not in w
        elif isinstance(dw,int):
            return dw > 0
        elif isinstance(w,float):
            print("obs: received w as float")
            return dw > 0

        return False


    try:

        if "scale" in json_data and "width" in json_data:
            w=json_data["width"]
            if valid_w(w):
                width=int(w)
                message+="Server received scale AND width, taking width."
            else:
                s=json_data["scale"]
                if valid_s(s):
                    scale=float(s)
                else:
                    missing.append(f"Server received scale AND width, both invalid:{s},{w}")

        elif "scale" in json_data and "width" not in json_data:
            s=json_data["scale"]
            if valid_s(s):
                scale=float(s)
            else:
                missing.append(f"The scale value received is not valid")
        elif "width" in json_data and "scale" not in json_data:
            w=json_data["width"]
            if valid_w(w):
                width=int(w)
                if width >MAX_IMAGE_WIDTH:
                    width=MAX_IMAGE_WIDTH
            else:
                missing.append(f"The width value received is not valid")
    except Exception as e:
        print("Exception determining scale of image for conversion: ",e)
        print(traceback.format_exc())
        return jsonify({"message":message,
            "error": "Error determining scale of image for conversion"}), 400

    if missing:
        return jsonify({"message": "".join(missing), "error": "Error during image conversion."}), 400

    try:
        #first to identify to get info on image levels
        responsedict={}
        stdout, stderr, error=identify_cmd(json_data)

        if error is not None:
            print("error from identify_cmd in convert was not none")
            raise error

        if stderr:
            stderrdec = stderr.decode()
            responsedict["stderr"] = stderrdec[:449] + "(...)" + stderrdec[-899:] if len(stderrdec) > 900 else stderrdec
            responsedict["message"] = "Error running identify command while trying to convert image."
            if not stdout:
                return jsonify(responsedict), 400

        stdec = stdout.decode()
        responsedict["stdout"] = stdec
        ld = parse_identify_str(stdec)
        responsedict["identify"]=json.dumps(ld)

        mainlevel = None
        wantedmax = MAX_IMAGE_WIDTH
        deal_w_scale_message = ""

        if "[unique]" in ld:
            mainlevel = "[unique]"
            if scale is not None:
                wantedmax = min(float(ld[mainlevel][0]) * scale, MAX_IMAGE_WIDTH)
                if wantedmax < MAX_IMAGE_WIDTH:
                    deal_w_scale_message += f"The scale proposed ({scale}) makes the final image larger than the maximum size. Defaulting to maximum size allowed."
            elif width is not None:
                wantedmax = min(width, MAX_IMAGE_WIDTH)
                if wantedmax < width:
                    deal_w_scale_message += f"The width proposed ({width}) is larger than the maximum size. Defaulting to maximum size allowed."
            elif width is None and scale is None:
                    wantedmax=MAX_IMAGE_WIDTH
                    deal_w_scale_message+= f"Width and scale received were None. Defaulting to maximum size allowed."
        else:
            mainlevel="[0]"
            if width is not None:
                for l in ld:
                    if ld[l][0] > width + 5:
                        mainlevel = l
                wantedmax = width
            elif scale is not None:
                wantedmax = min(ld["[0]"][0] * scale, MAX_IMAGE_WIDTH)
                if wantedmax < MAX_IMAGE_WIDTH:
                    deal_w_scale_message += f"The scale proposed ({scale}) makes the final image larger than the maximum size. Defaulting to maximum size allowed"
                for l in ld:
                    if ld[l][0] > wantedmax + 5:
                        mainlevel = l
            else:
                deal_w_scale_message += "Width and scale received were None. Defaulting to maximum size allowed"
                for l in ld:
                    if ld[l][0] > wantedmax + 5:
                        mainlevel = l

        print("sending this to convert_cmd:\n",json_data, ld, mainlevel, wantedmax)
        output, stderr, error=convert_cmd(json_data, ld, mainlevel, wantedmax)

        output["message"] = output.get("message", "") + message + deal_w_scale_message

        print(output)

        return jsonify(output), 200

    except Exception as e:
        print("Exception in convert: ",e)
        print(traceback.format_exc())
        return jsonify({"error": "Error running convert command"}), 400

def get_position_df(tp_path):
    position_df=None
    if os.path.isfile(tp_path):
        if "list" in tp_path:
            position_df = pd.read_csv(
                tp_path,
                names=["barcode", "in_tissue", "row", "col", "pxl_row_in_fullres", "pxl_col_in_fullres"],
            )
        else:
            position_df = pd.read_csv(
                tp_path,
                header=0,
                names=["barcode", "in_tissue", "row", "col",  "pxl_row_in_fullres", "pxl_col_in_fullres"],
            )

    position_df["in_tissue"]=position_df["in_tissue"].astype(bool)

    return position_df

def get_xy_barcode_dict(tp_path):
    position_df = get_position_df(tp_path)

    xy_keys = [
        f"{x}_{y}" for x, y in position_df[["pxl_col_in_fullres", "pxl_row_in_fullres"]].itertuples(index=False)
    ]

    barcode_values = position_df["barcode"].to_list()
    xy_barcode_dict = dict(zip(xy_keys, barcode_values))

    return xy_barcode_dict

def get_metric_summary_dict(ms_path):
    try:
        df = pd.read_csv(ms_path)
        metrics_dict = df.iloc[0].to_dict()
    except Exception as e:
        print("empty metric_summary_path", e)
        metrics_dict = {}
    return metrics_dict

def get_sample_metadata(metadata_path):
    scale_factors_dict = {"spot_diameter_fullres":300,"tissue_hires_scalef":0.5}

    if os.path.exists(metadata_path):
        with open(metadata_path, "r") as f:
            scale_factors_dict = json.load(f)
    else:
        print("Scale factors file is not available, resorting to defaults:",str(scale_factors_dict))

    return scale_factors_dict.get("spot_diameter_fullres", 0), scale_factors_dict.get(
        "tissue_hires_scalef", 0
    )


def create_gene_parquet(tp_path, fmh5_path, genes, is_ensembl):
    try:
        adata = sc.read_10x_h5(fmh5_path, gex_only=False)

        adata.var_names_make_unique()

        columns = adata.var["gene_ids"].to_list() if is_ensembl else adata.var.index

        columns = list(set(columns) & set(genes))

        if len(columns) < len(genes):
            print(
                f"Warning: The genes {set(genes) - set(columns)} were not found in the Spaceranger output and excluded from the parquet file"
            )

        if len(columns) == 0:
            print("WARNING: No genes were found to create the parquet file")

        # get tissue position info
        position_df = get_position_df(tp_path)
        position_df = position_df.set_index("barcode")

        common=list(set(adata.obs.index) & set(position_df.index.to_numpy()))
        if len(common)==0:
            raise Exception("The indices of the tissue positions file and the filtered_feature_bc_matrix.h5"+
                   "file are not the same, can't join them")

        print("sizes:",len(adata.obs.index),len(position_df.index.to_numpy()),len(common))

        # order adata observations to match that of positions
        adata = adata[common, :]
        adata.obs = adata.obs.join(position_df)

        # get x_y ids
        ids = [
            f"{x}_{y}" for x, y in adata.obs[[ "pxl_col_in_fullres","pxl_row_in_fullres"]].itertuples(index=False)
        ]

        if is_ensembl:
            df = pd.DataFrame(
                data=adata[:, adata.var["gene_ids"].isin(columns)].X.toarray(),
                index=ids,
                columns=columns,
            )
        else:
            df = pd.DataFrame(
                data=adata[:, columns].X.toarray(), index=ids, columns=columns
            )



        #I dont have a qc variable in the UI, if someone wants this it will have to be implemented by them
        # if qc:
        #     qc_df = sc.pp.calculate_qc_metrics(adata)[0]
        #     qc_df.index = ids
        #     df = pd.concat([df, qc_df], axis=1)

        df.index.name = "__index_level_0__"

        option=1 #or 1

        if option==0:

            print("attempting to save parquet then load as binary")

            uid=generate_random_string(10)
            parquet_file_path = f'/scratch/{uid}.parquet'
            df.to_parquet(parquet_file_path)

            binary_data=""

            # Step 2: Read the Parquet file in binary mode
            with open(parquet_file_path, 'rb') as file:
                binary_data = file.read()

            return binary_data

        else:
            binary_data = df.to_parquet()
            print("sending parquet as is")
            return binary_data


    except Exception as e:
        print("Error creating gene expression, exception as e: ",e)
        print(traceback.format_exc())
        return jsonify({"error": "Error creating gene expression"}), 500


@app.route('/upload_gene_data', methods=['POST'])
def upload_gene_data():
    if not request.is_json:
        return jsonify({"error": "Request is not json formatted"}), 400

    json_data=request.json

    error_messages={
        "token":"We need a token to authorize.",
        "rid":  "We need a role id.",
        "project":"We need a project name.",
        "sample": "We need a sample name.",
        "genes": "We need a list of genes.",
        "tissue_positions_loc": "We need the location of the tissue postiions file.",
        "filtered_matrix_h5_loc": "We need the location of the filtered matrix h5 file.",
        "is_ensembl": "We need to know if the genes have ensembl names."
    }

    missing=[error_messages[i] for i in error_messages if i not in json_data]

    #print(json_data)

    if len(missing)>0:
        request_json={}
        request_json["message"]="".join(missing)
        request_json["error"]="Missing input for gene upload"
        return jsonify( request_json ), 400
    #create_gene_parquet(tp_path, fmh5_path, genes, is_ensembl):
    binary_data = create_gene_parquet(json_data["tissue_positions_loc"],
                                      json_data["filtered_matrix_h5_loc"],
                                      json_data["genes"],
                                      json_data["is_ensembl"])

    token=json_data["token"]
    rid=json_data["rid"]
    project_id=json_data["project"]
    name=json_data["sample"]

    headers = {
        "Content-Type": "Content-Type: multipart/form-data",
        "Authorization": f"bearer {token}",
    }
    response = requests.post(
        CLOUD_ADDRESS + f"add_sample_parquet?rid={rid}&project_id={project_id}&sample={name}",
        headers=headers, data=binary_data
    )

    request_json= {"message": "From GCP about gene parquet: " + response.text}

    return jsonify( request_json ), response.status_code

def _add_distance_to_file_scores(possibilities):
    minl=np.inf;maxl=-np.inf
    for p in possibilities:
        if len(p["file"])<minl:
            minl=len(p["file"])
        if len(p["file"])>maxl:
            maxl=len(p["file"])

    dist=maxl-minl
    if dist==0:
        dist=1
    for p in possibilities:
        ind=(maxl-len(p["file"]))/dist
        if len(p["file"])==0:
            p["score"]=-1
        else:
            p["score"]+=ind

    return possibilities

def _return_best_fit(found, files):
    for n in found:
        #file to lookfor
        needle=n[0]
        #print("needle",needle)
        possibilities=[{"file":"","score":-1}]
        ready=False
        if "tissue_positions.csv" in needle:
            #there may be many files with this string, prefer those
            # that are shorter and have the words outs and spatial in them
            for p in [
                {"file": f,
                "score": ("outs" in f)+("spatial" in f)} for f in files if "tissue_positions.csv" in f]:
                possibilities.append(p)
                ready=True

        if ".gz" in needle:
            needle=needle.replace(".gz","")

        if not ready:
            for p in [{"file": f,
                       "score": 0} for f in files if needle in f]:
                possibilities.append(p)

        #for p in possibilities:
        #    print("    -",p)

        possibilities=_add_distance_to_file_scores(possibilities)
        max_score = max(possibilities, key=lambda x: x["score"])
        n[1]=max_score["file"]

    return found

@app.route('/find_10x_highres_image', methods=['POST'])
def find_10x_highres_image():
    if not request.is_json:
        return jsonify({"error": "Request is not json formatted"}), 400

    json_data=request.json

    error_messages={
        "sr_path":  "We need a spaceranger path.",
        "sample":"We need a sample name.",
        "scalefactors_json": "We need a scalefactors_json file"
    }

    missing=[error_messages[i] for i in error_messages if i not in json_data]

    if len(missing)>0:
        request_json={}
        request_json["message"]="".join(missing)
        request_json["error"]="Missing input for search"
        return jsonify( request_json ), 400

    received_scalefactors_loc=json_data["scalefactors_json"]
    print(f"received_scalefactors_loc {received_scalefactors_loc}")

    search_here=os.path.join(json_data["sr_path"],json_data["sample"])

    thesefiles=["tissue_hires_image.png"]

    if not os.path.exists(received_scalefactors_loc):
        thesefiles.append("scalefactors_json.json")

    found=[[n,""] for n in thesefiles]
    files = glob.glob(f'{search_here}/**/**', recursive = True)

    found=_return_best_fit(found,files)

    #image and scaleloc read
    finalimageloc=""
    finaljsonloc=received_scalefactors_loc
    finalscale=0

    for f in found:
        if "tissue_hires_image" in f[0] and len(f[1])>0:
            finalimageloc=f[1]
            continue
        if "scalefactors_json" in thesefiles and len(f[1])>0:
            finaljsonloc=f[1]
            continue
        else:
            finaljsonloc=received_scalefactors_loc


    message=""
    code=200

    scale_factors_dict = {}
    if len(finaljsonloc)>0:
        #load and get scale
        with open(finaljsonloc, "r") as f:
            scale_factors_dict = json.load(f)
        finalscale=scale_factors_dict.get( "tissue_hires_scalef", 0)
    else:
        message+="No scalefactors_json was found."

    if len(finalimageloc)==0:
        message+="No tissue_hires_image was found."
        if finalscale!=0:
            message+="However a scale was found inside the scalefactors_json.json file"
        code=400

    request_json={}

    request_json["message"]=message
    request_json["image_loc"]=finalimageloc
    request_json["image_scale"]=finalscale

    return jsonify( request_json ), code

@app.route('/find_srfiles', methods=['POST'])
def find_srfiles():
    if not request.is_json:
        return jsonify({"error": "Request is not json formatted"}), 400

    json_data=request.json

    error_messages={
        "sr_path":  "We need a spaceranger path.",
        "sample":"We need a sample name.",
        "file_list": "We need a file list to search for."
    }

    missing=[error_messages[i] for i in error_messages if i not in json_data]

    if len(missing)>0:
        request_json={}
        request_json["message"]="".join(missing)
        request_json["error"]="Missing input for search"
        return jsonify( request_json ), 400

    search_here=os.path.join(json_data["sr_path"],json_data["sample"])

    found=[[n,""] for n in json_data["file_list"]]
    files = glob.glob(f'{search_here}/**/**', recursive = True)

    #now it's going to be a special case for every file

    found=_return_best_fit(found,files)

    request_json={}

    request_json["message"]="searched for files correctly"
    request_json["file_locations"]=found

    return jsonify( request_json ), 200

@app.route('/initialize_in_db', methods=['POST'])
def initialize_in_db():
    if not request.is_json:
        return jsonify({"error": "Request is not json formatted"}), 400

    json_data=request.json
    cloud_address=CLOUD_ADDRESS

    if "cloud_address" in json_data:
        cloud_address=json_data["cloud_address"]

    error_messages={
        "token":  "We need a token to authorize.",
        "rid":  "We need a role id.",
        "project":"We need a project name.",
        "sample": "We need a sample name."
    }

    lookfor=["token","rid", "project", "sample" ]

    missing=[error_messages[i] for i in lookfor if i not in json_data]

    if len(missing)>0:
        request_json={}
        request_json["message"]="".join(missing)
        request_json["error"]="Missing input for search"
        return jsonify( request_json ), 400

    request_json = {
        "rid": json_data["rid"],
        "project_id": json_data["project"],
        "sample": json_data["sample"],
        "sample_fields": {},
        "xy_barcode": {},
    }

    token=json_data["token"]

    headers = {"Content-Type": "application/json", "Authorization": f"bearer {token}"}

    response = requests.post(
        cloud_address + "add_sample", headers=headers, json=request_json
    )

    request_json= {"message": "From GCP: " + response.text}

    return jsonify( request_json ), response.status_code



@app.route('/upload_spot_info', methods=['POST'])
def upload_spot_info():
    #if it gets to this point is because the image complies with the standards, try.
    if not request.is_json:
        return jsonify({"error": "Request is not json formatted"}), 400

    cloud_address="https://add-sample-barcodes-142858704207.us-central1.run.app/"
    error_messages={
        "token":  "Token is needed to authenticate request.",
        "project" :  "Project name is needed.",
        "sample":"We need a sample name",
        "rid"  :  "We need a role.",
        "tp_path": "A local path for tissue_positions.csv is needed."
    }

    warning_messages={
        "ms_path":  "A metrics summary path can be added to display metric information in the celery interface.",
        "metadata_path" :  "A scalefactors_json.json path can be usefull if uploading the 10X's 'higres_tissue image'"
    }

    json_data=request.json

    missing=[error_messages[i] for i in error_messages if i not in json_data]

    if len(missing)>0:
        request_json={}
        request_json["message"]="".join(missing)
        request_json["error"]="Missing input for gene upload"
        return jsonify( request_json ), 400

    token=json_data["token"]
    rid=json_data["rid"]
    project_id=json_data["project"]
    name=json_data["sample"]

    tp_path=json_data["tp_path"]

    to_warn=[warning_messages[i] for i in warning_messages if i not in json_data]

    metrics_summary_dict = {}
    spot_diameter_fullres=300
    tissue_hires_scalef=0
    sample_fields={}
    if "ms_path" in json_data:
        sample_fields["metrics_summary"]=get_metric_summary_dict(json_data["ms_path"])

    if "metadata_path" in json_data:
        spot_diameter_fullres, tissue_hires_scalef=get_sample_metadata(json_data["metadata_path"])
        sample_fields["spot_diameter_fullres"]= spot_diameter_fullres
        sample_fields[ "tissue_hires_scalef"]= tissue_hires_scalef

    xy_barcode_dict = get_xy_barcode_dict(tp_path)

    request_json = {
        "rid": json_data["rid"],
        "project_id": json_data["project"],
        "sample": json_data["sample"]
    }

    if len(sample_fields)>0:
        request_json["sample_fields"]=sample_fields

    print("request_json",request_json)

    request_json["xy_barcode"]=xy_barcode_dict

    headers = {"Content-Type": "application/json", "Authorization": f"bearer {token}"}

    response = requests.post(
        cloud_address, headers=headers, json=request_json
    )

    warn_mesage=""
    if len(to_warn)>0:
        warn_mesage="To note:"
        warn_mesage+="".join(warning_messages)
    request_json= {"message": "From GCP: " + response.text+warn_mesage}

    return jsonify( request_json ), response.status_code


@app.route('/upload_image', methods=['POST'])
def upload_image():
    #if it gets to this point is because the image complies with the standards, try.
    if not request.is_json:
        return jsonify({"error": "Request is not json formatted"}), 400

    json_data=request.json

    cloud_address=CLOUD_ADDRESS
    error_messages={
        "token":  "Token is needed to authenticate request.",
        "project" :  "Project name is needed.",
        "sample":"We need a sample name",
        "rid"  :  "We need a role.",
        "image_path": "A local path for the image is needed.",
        "scale": "A scale is needed even if it is 1.0"
    }

    missing=[error_messages[i] for i in error_messages if i not in json_data]

    if len(missing)>0:
        request_json={}
        request_json["message"]="".join(missing)
        request_json["error"]="Error during request for image upload"
        return jsonify( request_json ), 400

    token=json_data["token"]
    image_path=json_data["image_path"]
    rid=json_data["rid"]
    name=json_data["sample"]
    project_id=json_data["project"]
    scale=json_data["scale"]

    img_extension=image_path.split(os.sep)[-1].split(".")[-1]

    print(f"sending this to google add_sample_image?rid={rid}&project_id={project_id}&sample={name}&img_extension={img_extension}&scale_f={scale}" )

    binary_data = None
    with open(image_path, "rb") as f:
        binary_data = f.read()

    headers = {
        "Content-Type": "Content-Type: multipart/form-data",
        "Authorization": f"bearer {token}",
    }

    response = requests.post(
        cloud_address
        + f"add_sample_image?rid={rid}&project_id={project_id}&sample={name}&img_extension={img_extension}&scale_f={scale}",
        headers=headers,
        data=binary_data,
    )

    print(response.text)
    #print(dir(response))

    request_json={}

    request_json["message"]=response.text

    return jsonify( request_json ), response.status_code


@app.route('/identify', methods=['POST'])
def identify():
    if not request.is_json:
        return jsonify({"error": "Request is not json formatted"}), 400

    try:
        responsedict={}
        json_data=request.json
        stdout, stderr, error=identify_cmd(json_data)

        if error is not None:
            print("error from identify_cmd was not none")
            raise error

        if stderr:
            stderrdec=stderr.decode()
            if len(stderrdec)>900:
                responsedict["stderr"]=stderrdec[0:449]+"(...)"+stderrdec[-899:]
            else:
                responsedict["stderr"]=stderrdec
            responsedict["error"]="output of identify wrote to stderr, maybe warning?"
            print("stderr was not empty",responsedict["error"])
            responsedict["message"]="Error running identify command"
            if not stdout:
                return jsonify(responsedict), 400

        stdec=stdout.decode()
        responsedict["stdout"]=stdec
        ld=parse_identify_str(stdec)
        responsedict["identify"]=json.dumps(ld)

        return jsonify(responsedict), 200
    except Exception as e:
        print("Error processing identify, exception as e: ",e)
        print(traceback.format_exc())
        return jsonify({"error": "Error running identify command"}), 400

@app.route('/list_sr_folders', methods=['GET','POST'])
def list_sr_folders():
    if not request.is_json:
        return jsonify({"error": "Request is not json formatted"}), 400

    try:
        json_data=request.json
        srd=json_data["spaceranger_dir"]
        entries = os.listdir(srd)
        directories = [entry for entry in entries if os.path.isdir(os.path.join(srd, entry)) and entry[0] != "."]
        return jsonify({"directories": directories}), 200
    except Exception as e:
        print(e)
        print(traceback.format_exc())
        return jsonify({"error": "Error while getting sr directories"}), 500

@app.route('/add_project_do', methods=['POST'])
def add_project_do():

    cloud_address=CLOUD_ADDRESS
    error_messages={
        "token":  "Token is needed to authenticate request.",
        "name" :  "Project name is needed.",
        "rid"  :  "Project needs a role.",
        "classes":"Project needs Classes."
    }

    """
    For now, is json should be true because we are not sending any form multipart.
    In fact, it will never do, because we are not sending from this form, we are sending
    from the nygc server at the location given in this form
    """

    if not request.is_json:
        return jsonify({"error": "Request is not json formatted"}), 400

    json_data=request.json

    lookfor=["cloud_address", "token", "name", "rid", "emails", "classes", "extra"]
    found={i:json_data[i] for i in lookfor if i in json_data}

    if "cloud_address"  in found:
        cloud_address=found["cloud_address"]

    request_json={}
    missing=[]

    if "token" not in found:
        missing.append(error_messages["token"])

    if "name" not in found:
        missing.append(error_messages["name"])
    else:
        request_json["project_id"]=found["name"]

    if "rid" not in found:
        missing.append(error_messages["rid"])
    else:
        request_json["rid"]=found["rid"]

    if "classes" in found:
        if len(found["classes"])==0:
            missing.append(error_messages["classes"])
        else:
            request_json["annotations"]=found["classes"]
    else:
        missing.append(error_messages["classes"])

    if "emails" in found:
        if len(found["emails"])==0:
            missing.append(error_messages["emails"])
        else:
            request_json["user_emails"]=found["emails"]
    else:
        missing.append(error_messages["emails"])

    if "extra" in found:
        if len(found["extra"])>0:
            request_json["extra_fields"]={}
            for e in found["extra"]:
                key=e["name"]
                value=e["value"]
                request_json["extra_fields"][key]=value


    if len(missing)==0:
        token=found["token"]
        headers = {"Content-Type": "application/json", "Authorization": f"bearer {token}"}

        response = requests.post(
            cloud_address + "initialize_project", headers=headers, json=request_json
        )

        print(response.text)
        print(dir(response))

        request_json["message"]=response.text

        return jsonify( request_json ), response.status_code
    else:
        message="".join(missing)
        return jsonify({"error": message}), 400

@app.route('/add_image')
def contact():
    return render_template('contact.html')

if __name__ == '__main__':
    port=5000
    if len(sys.argv)>1:
        port=int(sys.argv[1])
    app.run(debug=True,host='0.0.0.0', port=port)