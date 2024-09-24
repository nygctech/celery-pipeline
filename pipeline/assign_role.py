import argparse
import requests
import config

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Assign a role to users' emails, even if the user has not signed-up yet."
    )
    parser.add_argument(
        help="rid (role id) to assign to users and initalize data bucket if it does not exist.",
        dest="rid",
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
    email_file = args.email_file
    token = args.token

    emails = []
    with open(email_file, "r") as f:
        emails = [l.strip() for l in f.readlines()]

    headers = {"Content-Type": "application/json", "Authorization": f"bearer {token}"}

    request_json = {"rid": rid, "user_emails": emails}

    response = requests.post(
        config.CLOUD_FN_URL + "new_role", headers=headers, json=request_json
    )
    print(response.text)
