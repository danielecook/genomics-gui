from flask import Flask, render_template, url_for
from flask_bootstrap import Bootstrap
from jinja2 import Template, environmentfilter
import os, sys

import yaml
import json
from tabulate import tabulate
import os
from yaml import load, dump
from pprint import pprint as pp
from collections import OrderedDict

import httplib2
import pprint

from apiclient.discovery import build
from apiclient.errors import HttpError
from apiclient.http import MediaIoBaseUpload
import io

from oauth2client.client import OAuth2WebServerFlow
from oauth2client.client import AccessTokenRefreshError
from oauth2client.client import flow_from_clientsecrets
from oauth2client.file import Storage
from oauth2client.tools import run

def connect():
    FLOW = flow_from_clientsecrets('client_secrets.json',
                                 scope=['https://www.googleapis.com/auth/bigquery',
                                        'https://www.googleapis.com/auth/genomics',
                                        'https://www.googleapis.com/auth/devstorage.read_write'])
    # Authorization and credentials
    storage = Storage('credentials.dat')
    credentials = storage.get()

    if credentials is None or credentials.invalid:
        credentials = run(FLOW, storage)

    http = httplib2.Http()
    http = credentials.authorize(http)

    bq = build('bigquery', 'v2', http=http)
    storage = build('storage','v1',http=http)
    genomics = build('genomics', 'v1beta2', http=http)
    return bq, storage, genomics

def submit(action):
  try:
    return action.execute()
  except HttpError as err:
    for error in json.loads(err.content)["error"]["errors"]:
      print error["reason"], error["message"]
  except AccessTokenRefreshError:
    print ("Credentials have been revoked or expired, please re-run"
           "the application to re-authorize")

app = Flask(__name__)
Bootstrap(app)

bq, storage, genomics = connect()
config = yaml.load(open("config.yaml",'r'))
PROJECT_NUMBER = config["PROJECT_NUMBER"]


def get_dataset_id_by_name(dataset_name):
    dataset_list = submit(genomics.datasets().list(projectNumber=PROJECT_NUMBER))["datasets"]
    dataset = [x for x in dataset_list if x["name"] == dataset_name]
    if len(dataset) == 1:
        return dataset[0]["id"]
    else:
        return None

@app.route("/")
def home():
    bq, storage, genomics = connect()
    datasets = submit(genomics.datasets().list(projectNumber=PROJECT_NUMBER))["datasets"]
    return render_template('view.html', **locals())

@app.route("/Project/")
def Project():
    datasets = submit(genomics.datasets().list(projectNumber=PROJECT_NUMBER))["datasets"]
    breadcrumbs = OrderedDict([("Project", "active")])
    return render_template('project.html', **locals())

@app.route("/Project/<dataset_name>/")
def Dataset(dataset_name):
    dataset_id = get_dataset_id_by_name(dataset_name)
    if dataset_id is None:
        raise Exception("Multiple databases with the same name not supported yet.")
    variantSets = submit(genomics.variantsets().search(body={"datasetIds":[dataset_id]}))
    variantSets = variantSets["variantSets"]
    # Reorganize variant sets
    metadata_display = ["fileformat", "fileDate", "source", "reference", "assembly" ]
    for n,varset in enumerate(variantSets):
        if "metadata" in varset:
            for meta in varset["metadata"]:
                for k,v in meta.items():
                    if v in metadata_display:
                        variantSets[n][meta["key"]] = meta["value"]
                        
    breadcrumbs = OrderedDict([("Project", url_for('Project')), (dataset_name, "active")])
    return render_template('dataset.html', **locals())

@app.route("/Project/<dataset_name>/<variant_set_id>/<ref>/")
def Variants(dataset_name, variant_set_id, ref):
    dataset_id = get_dataset_id_by_name(dataset_name)
    if dataset_id is None:
        raise Exception("Multiple databases with the same name not supported yet.")
    variantSets = submit(genomics.variantsets().search(body={"datasetIds":[dataset_id]}))
    variantSets = variantSets["variantSets"]
    breadcrumbs = OrderedDict([("Project", url_for('Project')),
                            (dataset_name, url_for('Dataset', dataset_name=dataset_name))])
    variants = submit(genomics.variants().search(body={"referenceName":ref, "variantSetIds":[variant_set_id]}))
    # Reorganize Callsets
    if "variants" in variants:
        sample_set = set()
        for k,v in enumerate(variants["variants"]):
            variants["variants"][k]["gt"] = {}
            for c in v["calls"]:
                sample_set.add(c["callSetName"])
                sample_name = c["callSetName"]
                gt_map = {k+1:v for k,v in enumerate(v["alternateBases"])}
                gt_map[0] = v["referenceBases"]
                variants["variants"][k]["gt"][sample_name] = c["genotype"]
        sample_set = sorted(sample_set)

    return render_template('variant.html', **locals())

@app.route("/<chrom>/<start>/<end>/")
def query(chrom, start, end):
    variants = vcf.Reader(filename="ANNOTATIONS.bcf").fetch(chrom,start,end)
    return render_template('index.html', **locals())


if __name__ == "__main__":
    app.run(debug=True)
    # Enter your Google Developer Project number
