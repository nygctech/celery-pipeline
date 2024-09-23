import argparse
import requests
import os
import pathlib
import pandas as pd
import sys
import scanpy as sc
import io
import warnings
import json
from io import BytesIO
import PIL
import math
import config

def get_position_df(spaceranger_path):
    """
    Reads outs/spatial/tissue_positions[_list].csv from the spaceranger output
    and returns the DataFrame with in_tissue barcodes.
    """
    if os.path.isfile(os.path.join(spaceranger_path, "outs/spatial/tissue_positions_list.csv")):
        position_df = pd.read_csv(os.path.join(spaceranger_path, "outs/spatial/tissue_positions_list.csv"), 
                                  names=['barcode', 'in_tissue', 'row', 'col', 'imageY', 'imageX'])
    else:
        position_df = pd.read_csv(os.path.join(spaceranger_path, "outs/spatial/tissue_positions.csv"), header=0,
                                  names=['barcode', 'in_tissue', 'row', 'col', 'imageY', 'imageX'])
    
    position_df = position_df[position_df['in_tissue'].astype(bool)]

    return position_df

def get_xy_barcode_dict(spaceranger_path):
    """
    Reads outs/spatial/tissue_positions[_list].csv from the spaceranger output.

    Returns a dict of x_y coordinate of spots as keys and the spot barcodes as values.
    """
    position_df = get_position_df(spaceranger_path)

    xy_keys = [f"{x}_{y}" for x, y in position_df[['imageX', 'imageY']].itertuples(index=False)]

    barcode_values = position_df['barcode'].to_list()
    xy_barcode_dict = dict(zip(xy_keys, barcode_values))
    
    return xy_barcode_dict

def get_metric_summary_dict(spaceranger_path):
    """
    Reads the two row metrics_summary.csv file and return it as a dict
    """
    try:
        metrics_path = os.path.join(spaceranger_path, "outs/metrics_summary.csv")
        df = pd.read_csv(metrics_path)
        metrics_dict = df.iloc[0].to_dict()
    except FileNotFoundError as e:
        print(e, sys.stderr)
        metrics_dict = {}
    return metrics_dict

def handle_server_response(response):
    """
    Prints the text in the response and raises exception if status code isn't 2xx
    """
    print(f"[Server response] {response.text}")
    response.raise_for_status()

def get_sample_metadata(sample_path):
    """
    Get spot_diameter_fullres, tissue_hires_scalef from sample
    """
    scale_factors_path = os.path.join(sample_path, "outs/spatial/scalefactors_json.json")
    scale_factors_dict = {}

    with open(scale_factors_path, 'r') as f:
        scale_factors_dict = json.load(f)

    return scale_factors_dict.get('spot_diameter_fullres', 0), scale_factors_dict.get('tissue_hires_scalef', 0)


def initialize_sample(rid, project_id, sample_path, token, name):
    """
    Create the document for the sample and uploads the xy_barcode json to storage
    """
    metrics_summary_dict = get_metric_summary_dict(sample_path)
    xy_barcode_dict = get_xy_barcode_dict(sample_path)
    spot_diameter_fullres, tissue_hires_scalef = get_sample_metadata(sample_path)

    request_json = {
        "rid": rid,
        "project_id": project_id,
        "sample": name,
        "sample_fields": {
            'metrics_summary': metrics_summary_dict,
            'spot_diameter_fullres': spot_diameter_fullres,
            'tissue_hires_scalef': tissue_hires_scalef
        },
        "xy_barcode": xy_barcode_dict
    }

    headers = {
        "Content-Type": "application/json",
        "Authorization": f"bearer {token}"
    }

    response = requests.post(config.CLOUD_FN_URL + "add_sample", headers=headers, json=request_json)
    handle_server_response(response)



def add_sample_image(rid, project_id, sample_path, token, name):
    """
    Uploads the hi res image from the spaceranger output to the storage bucket
    """
    img_path = os.path.join(sample_path, 'outs/spatial/tissue_hires_image.png')
    img_extension = 'png'

    binary_data = None
    with open(img_path, 'rb') as f:
        binary_data = f.read()

    headers = {
        "Content-Type": "Content-Type: multipart/form-data",
        "Authorization": f"bearer {token}"
    }


    response = requests.post(config.CLOUD_FN_URL + f"add_sample_image?rid={rid}&project_id={project_id}&sample={name}&img_extension={img_extension}",
                            headers=headers, data=binary_data)
    handle_server_response(response)

def add_resized_fullres(rid, project_id, fullres_path, token, name):
    """
    Add a resized version of the fullres img which is smaller than config.MAX_IMG_SIZE and add it to the sample's document
    """
    PIL.Image.MAX_IMAGE_PIXELS = None
    img = PIL.Image.open(fullres_path)
    img_size = os.path.getsize(fullres_path)

    width, height = img.size   
    init_size_scale = math.sqrt(config.MAX_IMG_SIZE / img_size)

    crop_width = round(width * init_size_scale)
    crop_height = round(height * init_size_scale)

    crop_img_size = img_size

    # resize at least once (event if smaller than max size)
    # decrease width and height by a factor of 2/3 until the JPEG is under config.MAX_IMG_SIZE
    while crop_img_size > config.MAX_IMG_SIZE or crop_img_size == img_size:
        # load bytes into buffer
        crop_img_bytes = BytesIO()
        crop_img = img.resize((crop_width, crop_height))
        crop_img.save(crop_img_bytes, 'JPEG', quality=95)
        crop_img_size = crop_img_bytes.tell()

        init_size_scale *= 2/3
        crop_width = round(crop_width * init_size_scale)
        crop_height = round(crop_height * init_size_scale)
        crop_img_bytes.seek(0)

    # read and close buffer
    binary_data = crop_img_bytes.read()
    crop_img_bytes.close()
    scale_f = round(crop_img.size[0] / width, 8)
    headers = {
                "Content-Type": "Content-Type: multipart/form-data",
                "Authorization": f"bearer {token}"
    }

    img_extension = 'JPEG'
    response = requests.post(config.CLOUD_FN_URL + f"add_sample_image?rid={rid}&project_id={project_id}&sample={name}&img_extension={img_extension}&scale_f={scale_f}",
                            headers=headers, data=binary_data)
    handle_server_response(response)

def add_sample_patches(rid, project_id, fullres_path, token, name):
    """Creates multiple patches of fullres image (to not exceed 32MB) and uploads them"""
    PIL.Image.MAX_IMAGE_PIXELS = None

    fullres_size = os.path.getsize(fullres_path)  
    fullres_extension = pathlib.Path(fullres_path).suffix.strip('.')

    assert fullres_extension.upper() in ['PNG', 'JPG', 'JPEG'], 'The fullres image extension must be .PNG, .JPEG, or .JPG'
    pil_type = fullres_extension
    if pil_type.upper() == 'JPG':
        pil_type = 'JPEG'

    #approximate num of patches by dividing num of bytes by 20MB (and rounding up)
    num_patches = fullres_size // 20_000_000 + 1

    with PIL.Image.open(fullres_path) as fullres_image:
        width, height = fullres_image.size

        patch_width_decimal = width / num_patches
        

        for i in range(num_patches):
            left = i * math.floor(patch_width_decimal)
            upper = 0
            if i < num_patches - 1:
                right = i * math.floor(patch_width_decimal) + math.floor(patch_width_decimal) 
            else:
                right = width
            lower = height

            patch = fullres_image.crop((left, upper, right, lower))
            
            # start quality at 95, decrease until image size gets under PATH_MAX_SIZE
            patch_size = math.inf
            quality = 100 
            while patch_size > config.PATCH_MAX_SIZE:
                # highest quality starts at 95 (for jpeg)
                quality -= 5
                im_bytes = BytesIO()
                patch.save(im_bytes, pil_type, quality=quality)
                patch_size = im_bytes.tell()
                im_bytes.seek(0)


            assert quality > 50, "Was unable to get fullres image patches to high enough quality, try scaling down the fullres image"
            
    
            binary_data = im_bytes.read()
            im_bytes.close()

            headers = {
                "Content-Type": "Content-Type: multipart/form-data",
                "Authorization": f"bearer {token}"
            }

            response = requests.post(config.CLOUD_FN_URL + f"add_sample_patch?rid={rid}&project_id={project_id}&sample={name}&img_extension={fullres_extension}&current_patch={i+1}&num_patches={num_patches}",
                                    headers=headers, data=binary_data)
            handle_server_response(response)


def add_sample_parquet(rid, project_id, sample_path, token, genes, is_ensembl, name, qc=False):
    """
    Creates a parquet file to upload to the storage, only including the specified genes
    """
    outs_path = os.path.join(sample_path, "outs")

    #reload the adata object as a 10x_h5 to grab features that aren't genes after getting ids
    if os.path.exists(os.path.join(outs_path, "custom_features.h5ad")):
        adata = sc.read_h5ad(os.path.join(outs_path, "custom_features.h5ad"))
    else:
        adata = sc.read_10x_h5(os.path.join(outs_path, "filtered_feature_bc_matrix.h5"), gex_only=False)

    adata.var_names_make_unique()

    columns = adata.var['gene_ids'].to_list() if is_ensembl else adata.var.index

    columns = list(set(columns) & set(genes))
    
    if len(columns) < len(genes):
        print(f"Warning: The genes {set(genes) - set(columns)} were not found in the Spaceranger output and excluded from the parquet file")


    if len(columns) == 0:
        print("WARNING: No genes were found to create the parquet file")

    #get tissue position info
    position_df = get_position_df(sample_path)
    position_df = position_df.set_index('barcode')

    #order adata observations to match that of positions
    adata = adata[position_df.index, :]
    adata.obs = adata.obs.join(position_df)

    #get x_y ids
    ids = [f"{x}_{y}" for x, y in adata.obs[['imageX', 'imageY']].itertuples(index=False)]

    if is_ensembl:
        df = pd.DataFrame(data = adata[:, adata.var['gene_ids'].isin(columns)].X.toarray(), index=ids, columns=columns)
    else:
        df = pd.DataFrame(data = adata[:, columns].X.toarray(), index=ids, columns=columns)

    if qc:
        qc_df = sc.pp.calculate_qc_metrics(adata)[0]
        qc_df.index = ids
        df = pd.concat([df, qc_df], axis=1)

    df.index.name = '__index_level_0__'
    binary_data = df.to_parquet()

    headers = {
        "Content-Type": "Content-Type: multipart/form-data",
        "Authorization": f"bearer {token}"
    }

    response = requests.post(config.CLOUD_FN_URL + f"add_sample_parquet?rid={rid}&project_id={project_id}&sample={name}",
                            headers=headers, data=binary_data)
    handle_server_response(response)


def add_sample_annotation(rid, project_id, name, token, annotation_path):
    headers = {
        "Content-Type": "Content-Type: multipart/form-data",
        "Authorization": f"bearer {token}"
    }
    
    df = pd.read_csv(annotation_path)
    df = df.dropna()
    binary_data = df.to_csv(index=False)

    response = requests.post(config.CLOUD_FN_URL + f"add_sample_annotation?rid={rid}&project_id={project_id}&sample={name}",
                            headers=headers, data=binary_data)
    
    handle_server_response(response)


def add_sample(rid, project_id, sample_path, token, genes, is_ensembl, name=None, qc=False, fullres_path=None, resize_fullres=False, upload_fullres=False, annotation_path=None):
    #default name of sample to name of last directory in `sample_path`
    if name is None:
        name = pathlib.PurePath(sample_path).name
    
    warnings.simplefilter(action='ignore', category=UserWarning)
    print(f"Adding sample from {sample_path} with name {name}...")
    initialize_sample(rid, project_id, sample_path, token, name)
    if upload_fullres:
        print(f"Uploading fullres image from {fullres_path} as multiple patches (backend will stitch them back together)...")
        add_sample_patches(rid, project_id, fullres_path, token, name)

    if resize_fullres:
        print(f"Resizing the fullres image to upload, this may take a while...")
        add_resized_fullres(rid, project_id, fullres_path, token, name)
    else:
        print(f"Uploading image from {sample_path} with name {name}...")
        add_sample_image(rid, project_id, sample_path, token, name)

    if len(genes) == 0:
        print("WARNING: no genes were listed so no Parquet file will be uploaded")
    else:
        print(f"Processing and uploading parquet file (gene expression) from {sample_path}...")
        add_sample_parquet(rid, project_id, sample_path, token, genes, is_ensembl, name, qc=qc)

    if annotation_path is not None:
        print(f"Uploading annotations from {annotation_path}...")
        add_sample_annotation(rid, project_id, name, token, annotation_path)

    warnings.simplefilter(action='default')

def main():
    parser = argparse.ArgumentParser(description='Create a project and process the samples')
    parser.add_argument(help='rid (role id) to assign to users and initalize data bucket if it does not exist.', dest='rid')
    parser.add_argument(help='Name of project.', dest='project_id')
    parser.add_argument(help='Paths to Spaceranger objects', dest='sample_paths', nargs='+')
    parser.add_argument(help='Path to file with names of genes on each line which will be used to to construct gene expression matrix', dest='gene_path')
    parser.add_argument('-n', '--names', help='Optional alternative names of samples - default names are the Spaceranger subdirectory names', dest='sample_names', nargs='+')
    parser.add_argument('-e', '--ensembl', help='Names of genes in gene_path are ensembl codes', dest='is_ensembl', action='store_true')
    parser.add_argument('-q', '--qc_metrics', help='Calculate spot QC metrics for the samples and add to the expression profile.', dest='qc', action='store_true')
    parser.add_argument('-t', '--token', help='token of current user to use for authentication.', dest='token', required=True)
    parser.add_argument('-f', '--fullres_paths', help='Path to the fullres images.', dest='fullres_paths', nargs='+')
    parser.add_argument('-u', '--upload_fullres', help='Upload the fullres image.', dest='upload_fullres', action='store_true')
    parser.add_argument('-c', '--resize_fullres', help=f'Upload a resized version of the fullres image which doesn\'t exceed {config.MAX_IMG_SIZE}.', dest='resize_fullres', action='store_true')
    parser.add_argument('-a', '--annotation_paths', help='Path to existing annotations', dest='annotation_paths', nargs='+')

    args = parser.parse_args()
    
    rid = args.rid
    project_id = args.project_id
    sample_paths = args.sample_paths
    sample_names = [None] * len(sample_paths) if args.sample_names is None else args.sample_names
    token = args.token
    gene_path = args.gene_path
    is_ensembl = args.is_ensembl
    qc = args.qc
    fullres_paths = [None] * len(sample_paths) if args.fullres_paths is None else args.fullres_paths #default list of None
    resize_fullres = args.resize_fullres
    upload_fullres = args.upload_fullres
    annotation_paths = [None] * len(sample_paths) if args.annotation_paths is None else args.annotation_paths #default list of None



    assert sample_names is None or len(sample_names) == len(sample_paths), f"The number of arguments in --names should be the the same as the number of Spaceranger paths specified ({len(sample_paths)})"
    assert fullres_paths is None or len(fullres_paths) == len(sample_paths), f"The number of arguments in --fullres_paths should be the the same as the number of Spaceranger paths specified ({len(sample_paths)})"
    assert annotation_paths is None or len(annotation_paths) == len(sample_paths), f"The number of arguments in --annotations_paths should be the the same as the number of Spaceranger paths specified ({len(sample_paths)})"

    if resize_fullres or upload_fullres:
        assert None not in fullres_paths, f"fullres_paths must be provided if the resize_fullres or upload_fullres flags are present"

    genes = []
    with open(gene_path, 'r') as f:
        genes = [l.strip() for l in f.readlines()]

    for sample_path, sample_name, fullres_path, annotation_path in zip(sample_paths, sample_names, fullres_paths, annotation_paths):
        add_sample(rid, project_id, sample_path, token, genes, is_ensembl, name=sample_name, qc=qc, fullres_path=fullres_path,
                    resize_fullres=resize_fullres, upload_fullres=upload_fullres, annotation_path=annotation_path)

if __name__ == '__main__':
    main()
