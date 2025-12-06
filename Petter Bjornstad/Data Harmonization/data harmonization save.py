# Libraries
import os
import sys
import pandas as pd
import boto3
import json
import getpass
import tempfile
import shutil
import hashlib
import base64
from botocore.client import Config
from boto3.s3.transfer import TransferConfig

user = getpass.getuser() 

sys.path.insert(0, os.path.expanduser('~') +
                  "/GitHub/CHCO-Code/Petter Bjornstad/Data Harmonization")
from data_harmonization import harmonize_data

if user == "choiyej":
    base_data_path = "/Users/choiyej/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Bjornstad/Biostatistics Core Shared Drive/"
    keys_path = ""
elif user == "pylell":
    base_data_path = "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive/"
    keys_path = ""
elif user == "shivaniramesh":
    base_data_path = os.path.expanduser("~/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/")
    keys_path = "/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/keys.json"
else:
    sys.exit(f"Unknown user: please specify root path for this user. (Detected user: {user})")

with open(keys_path, "r") as f:
    keys = json.load(f)
session = boto3.Session(
    aws_access_key_id=keys['MY_ACCESS_KEY'],
    aws_secret_access_key=keys['MY_SECRET_KEY'],
)
s3 = session.client(
    "s3",
    endpoint_url="https://s3.kopah.uw.edu",
    config=Config(
        signature_version="s3v4",
        s3={"use_accelerate_endpoint": False},
        retries={"max_attempts": 3},
    )
)

clean = harmonize_data()
clean.to_csv(base_data_path + "Data Harmonization/Data Clean/harmonized_dataset.csv", index=False)


tmp_path = "/tmp/harmonized_upload.csv"
clean.to_csv(tmp_path, index=False, lineterminator="\n", encoding="utf-8")

print("Temporary CSV path:", tmp_path)
print("File size:", os.path.getsize(tmp_path))

# Upload
try:
    s3.upload_file(
        Filename=tmp_path,
        Bucket="harmonized.dataset",
        Key="harmonized_dataset.csv",
    )
    print("Upload successful!")
except Exception as e:
    print(f"Upload failed: {e}")
finally:
    if os.path.exists(tmp_path):
        os.remove(tmp_path)
        print("Temporary file removed.")