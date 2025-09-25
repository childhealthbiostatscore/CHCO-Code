# from https://www.notion.so/Setting-up-Kopah-1933695924ef80b187d7c2df978ece92

library(reticulate)
reticulate::use_python("/mmfs1/gscratch/scrubbed/mdanto/pytorch-cuda11/bin/python") # replace with your UWNetID
library(jsonlite) 

## Load boto3 and pandas
boto3 <- reticulate::import("boto3")
pd <- reticulate::import("pandas")

## Create an S3 client
keys <- fromJSON("/mmfs1/home/mdanto/keys.json") # replace with your UW ID
session <- boto3$session$Session(
  aws_access_key_id = keys$MY_ACCESS_KEY,
  aws_secret_access_key = keys$MY_SECRET_KEY
)

## Create an S3 client with the session
s3 <- session$client("s3", endpoint_url = "https://s3.kopah.uw.edu")

bucket = "today" # bucket name in Kopah
