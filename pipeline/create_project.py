import argparse
import requests
import config

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Creates a new Visium project")
    parser.add_argument(
        help="rid (role id) to assign to users and initalize data bucket if it does not exist.",
        dest="rid",
    )
    parser.add_argument(help="Name of project.", dest="project_id")
    parser.add_argument(
        help="Path to file with names of annotations on each line",
        dest="annotation_file",
    )
    parser.add_argument(help="Path to file with emails on each line", dest="email_file")
    parser.add_argument(
        "-t",
        "--token",
        help="token of current user to use for authentication.",
        dest="token",
        required=True,
    )

    args = parser.parse_args()
    rid = args.rid
    project_id = args.project_id
    annotation_file = args.annotation_file
    email_file = args.email_file
    token = args.token

    annotations = []
    with open(annotation_file, "r") as f:
        annotations = [l.strip() for l in f.readlines()]
    annotations = list(dict.fromkeys(annotations))  # remove duplicates, keep order

    emails = []
    with open(email_file, "r") as f:
        emails = [l.strip() for l in f.readlines()]

    headers = {"Content-Type": "application/json", "Authorization": f"bearer {token}"}

    request_json = {
        "rid": rid,
        "project_id": project_id,
        "annotations": annotations,
        "user_emails": emails,
        "extra_fields": {"isVisium": True},
    }

    response = requests.post(
        config.CLOUD_FN_URL + "initialize_project", headers=headers, json=request_json
    )
    print(response.text)
