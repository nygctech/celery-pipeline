import os

samples=["CU107-U54-HRA-101-A", "CU107-U54-HRA-101-B", "CU107-U54-HRA-113-C", "CU107-U54-HRA-113-D", "CU108-U54-HRA-101-A", "CU108-U54-HRA-101-B", "CU108-U54-HRA-113-C", "CU108-U54-HRA-113-D", "CU117-U54-HRA-127-C", "CU117-U54-HRA-127-D", "CU118-U54-HRA-127-C", "CU118-U54-HRA-127-D", "CU118-U54-HRA-182-A", "CU118-U54-HRA-182-B", "CU121-U54-HRA-185-C", "CU121-U54-HRA-185-D", "CU121-U54-HRA-192-A", "CU121-U54-HRA-192-B", "CU122-U54-HRA-185-C", "CU122-U54-HRA-185-D", "CU122-U54-HRA-192-A", "CU122-U54-HRA-192-B", "CU123-U54-HRA-097-C", "CU123-U54-HRA-097-D", "CU123-U54-HRA-130-A", "CU123-U54-HRA-130-B", "CU124-U54-HRA-097-C", "CU124-U54-HRA-097-D", "CU124-U54-HRA-130-A", "CU124-U54-HRA-130-B", "CU141-U54-HRA-216-C", "CU141-U54-HRA-216-D", "CU142-U54-HRA-216-C", "CU142-U54-HRA-216-D", "CU167-U54-HRA-284-A", "CU167-U54-HRA-284-B", "CU168-U54-HRA-284-A", "CU168-U54-HRA-284-B", "CU169-U54-HRA-290-C", "CU169-U54-HRA-290-D", "CU169-U54-HRA-296-A", "CU169-U54-HRA-296-B", "CU170-U54-HRA-290-C", "CU170-U54-HRA-290-D", "CU170-U54-HRA-296-A", "CU170-U54-HRA-296-B", "CU177-U54-HRA-311-C", "CU177-U54-HRA-311-D", "CU177-U54-HRA-317-A", "CU177-U54-HRA-317-B", "CU178-U54-HRA-311-D", "CU178-U54-HRA-317-A", "CU179-U54-HRA-229-A", "CU179-U54-HRA-229-B", "CU179-U54-HRA-305-C", "CU179-U54-HRA-305-D", "CU180-U54-HRA-229-A", "CU180-U54-HRA-229-B", "CU180-U54-HRA-305-C", "CU180-U54-HRA-305-D", "CU189-U54-HRA-317-C", "CU189-U54-HRA-317-D", "CU190-U54-HRA-317-C", "CU190-U54-HRA-317-D"]

imloc="/gpfs/commons/groups/phatnani_lab/mxia/u54_lsc_fullres_celery/"
project="U54_LSC"

samplenames={}
for s in samples:
    p1=s.split("-")[0]
    pl=s.split("-")[-1]
    samplenames[p1+pl]=s

imagenames={}
for f in os.listdir(imloc):
    p1=f.split("_")[0]
    pl=f.split("-")[-1][:-5]
    imagenames[p1+pl]=f


for s in samplenames:
    try:
        n=samplenames[s]
        im=imagenames[s]
        astr=""
        astr+=f"gcloud storage cp {imloc}{im} "
        astr+=f"gs://techinno.appspot.com/CGND/"
        astr+=f"{project}/{n}/{n}_fullres_image.jpg"
        print(astr)
    except:
        print(f"#error {s}")

srloc="/gpfs/commons/groups/phatnani_lab/Spatial_Multiomics/spaceranger_output/raw_spaceranger_output/"

for s in samplenames:
    try:
        n=samplenames[s]
        im=imagenames[s]
        astr=""
        astr+=f"gcloud storage cp {imloc}{im} "
        astr+=f"gs://techinno.appspot.com/CGND/"
        astr+=f"{project}/{n}/{n}_fullres_image.jpg"
        print(astr)
        print(f"gcloud storage cp {srloc}{n}/outs/filtered_feature_bc_matrix/* gs://techinno.appspot.com/CGND/{project}/{n}/filtered_feature_bc_matrix")
    except:
        print(f"#error {s}")

for s in samplenames:
    try:
        n=samplenames[s]
        im=imagenames[s]
        print(f"gcloud storage cp {srloc}{n}/outs/spatial/tissue_positions.csv gs://techinno.appspot.com/CGND/{project}/{n}/spatial/tissue_positions.csv")
        print(f"gcloud storage cp {srloc}{n}/outs/spatial/scalefactors_json.json gs://techinno.appspot.com/CGND/{project}/{n}/spatial/scalefactors_json.json")
    except:
        print(f"#error {s}")