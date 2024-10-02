# Description

This pipeline is designed to upload data from data sources (i.e. local machines, clusters, etc) to projects within the [Celery](https://techinno.firebaseapp.com/) webapp. 

This repository contains Python scripts which can be downloaded (see [Installation](#Installation)) and executed. Each script requires an additional `-t` flag followed by the user's current authentication token. This token can be copied to the user's clipboard by going to the Celery website (https://techinno.firebaseapp.com/) and pressing the "Get Token" button at the bottom.

# Installation

Clone the repository into the current directory:

```
$ git clone https://github.com/nygctech/celery-pipeline.git
```

Install the required Python packages from an environment with Python>=3.9:

```
$ pip install -r celery-pipeline/requirements.txt
```

# Adding new users

Roles dictate which projects users may be added to, each project is associated with a particular role, and multiple projects may be associated with that role. To access a specific project however, a user must still be added to that project (see [Create Visium project](#create-visium-project))

A user may be added with a role from an existing user with that role (see [below](#assigning-roles-to-new-and-existing-users)), or an admin may grant a new role to the user.

Notes on roles:
- Roles are can be thought more so as "groups" as opposed to having certain permissions like admins vs users
- Users can have multiple roles associated with their account
- In some of the outputs of the scripts, roles may be referred to as "rid" (role ID)


Currently, all users must use their @nygenome.org account to sign-in (this may be changed by an admin in the processSignUp Cloud Function).

**Known bug:** After a user signs-up for the first time in Celery, they may need to wait a few minutes before their permissions (specifically roles) get propagated and updated correctly.  

## Assigning roles to new and existing users

A user that has already been assigned a role may grant that same role to existing or new users. This can be done by first creating a text file with the users' emails on each line (named `emails.txt` for this example). Then, run the following command from the directory containing `celery-pipeline`:

```
$ python pipeline/assign_role.py [role name] emails.txt -t [token from Celery]
```

# Using the pipeline with Visium projects

## Create Visium project

A user can create a new Visium project, assign other users to the project, and specify the annotation region names for the project. This can be done by first creating a text file with the names of the annotations that will be used on each line (e.g. `annotations.txt`), and then creating a text file with the emails of users that will be performing the annotations (e.g. `projectUsers.txt`).

Then, run the following command from the directory containing `celery-pipeline`:

```
$ python pipeline/create_project [role name] [project name] annotations.txt projectUsers.txt -t [token from Celery]
```

A user may add or delete annotation names by revising the annotations text file, and re-running the above script. However, spots with annotation names that are deleted should also be deleted from the Celery viewer. Similarly, you may add or remove users from the project by revising the project users text file and re-running the script.

_Warning_ - Don't use an existing project name. If this does happen, the annotation names and users may get overwritten, so redo the annotation names and the users and re-run the script.

## Add Visium sample

A user may either upload a single sample ran through Spaceranger one at a time, or they may upload multiple samples at once. In both cases, a user must first create a text file which lists the names of marker genes on each line that they will use for annotation. 

Then, to upload a sample, they must run the following command:

```
$ python pipeline/add_sample.py [role name] [project name] [path to sample's Spaceranger folder] [path to gene txt file] \
    -f [path to image used for alignment] \
    --resize_fullres \
    --qc_metrics \
    -t [token from Celery]
```

The `-f [...]` and `--resize_fullres` flags may be omitted if the "hires" image generated by Spaceranger has enough detail to use for annotation. 

The full list of flags are their descriptions are:

```
positional arguments:
  rid                   rid (role id) to assign to users and initalize data bucket if it does not exist.
  project_id            Name of project.
  sample_paths          Paths to Spaceranger objects
  gene_path             Path to file with names of genes on each line which will be used to to construct gene expression matrix

optional arguments:
  -h, --help            show this help message and exit
  -n SAMPLE_NAMES [SAMPLE_NAMES ...], --names SAMPLE_NAMES [SAMPLE_NAMES ...]
                        Optional alternative names of samples - default names are the Spaceranger subdirectory names
  -e, --ensembl         Names of genes in gene_path are ensembl codes
  -q, --qc_metrics      Calculate spot QC metrics for the samples and add to the expression profile.
  -f FULLRES_PATHS [FULLRES_PATHS ...], --fullres_paths FULLRES_PATHS [FULLRES_PATHS ...]
                        Path to the fullres images.
  -u, --upload_fullres  Upload the fullres image.
  -c, --resize_fullres  Upload a resized version of the fullres image which doesn't exceed 30000000.
  -a ANNOTATION_PATHS [ANNOTATION_PATHS ...], --annotation_paths ANNOTATION_PATHS [ANNOTATION_PATHS ...]
                        Path to existing annotations
  -t TOKEN, --token TOKEN
                        token of current user to use for authentication.
```

# Other project types

An admin may need to initialize IF and HD projects. See the readme in the [Celery](https://github.com/mssanjavickovic/celery/) repository for further information on how to structure these projects.

# HD Visium

HD Visium projcets may be annotated in a similar method to other "HD" Celery projects. The fullres microscope image may be uploaded to an HD project (by an admin), and the user can created polygon-style annotations around regions. For best results, the microscope image should be cropped as much as possible before being aligned with Loupe Browser, which will make the image more manageable in Celery (the image may also be downsized again before being uploaded so that it doesn't take too long to load in the browser). 

The annotations can be exported from Celery into a `.json` file, and converted to a csv with barcodes in one columns and AARs in another column. To do this, run:

```
$ python pipeline/visium_hd_conversion.py [sample spaceranger path] [json annotations] [target bin size] -t [token from Celery]
```

The full list of flags are their descriptions are:

```
positional arguments:
  spaceranger_path      Path to the Spaceranger output files
  annotation_path       Path to the json annotation file exported by Celery
  target_bin_size       Target bin size, will automatically rebin if size does not exist

optional arguments:
  -h, --help            show this help message and exit
  --out_path OUT_PATH, -o OUT_PATH
                        Path to csv file with barcodes and AARs (default "AARs.csv")
  --save_parquet, -s    Save the rebinned data as a parquet file
  --xy_translation x y, -t x y
                        Translation applied to annotations in x and y directions (default 0 0)
  --scalefactor SCALEFACTOR, -f SCALEFACTOR
                        Scale factor applied to annotations after translation (default 1)
  --distance_threshold DISTANCE_THRESHOLD, -d DISTANCE_THRESHOLD
                        Number of barcodes away (in target bin size) a blank spot must be from an annotated region to
                        automatically carry-over nearby annotations (default 32).
```