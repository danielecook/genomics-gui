from flask import Flask, render_template, url_for, request, jsonify
from flask_bootstrap import Bootstrap
from jinja2 import Template, environmentfilter
import os, sys
import time, datetime

import yaml
import json
from json import loads
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

from variables import metadata_keys
from pprint import pprint as pp

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



@app.template_filter('fixtime')
def _jinja2_filter_datetime(timestamp, fmt=None):
    return datetime.datetime.fromtimestamp(float(timestamp)/1000).strftime("%Y-%m-%d %H:%M:%S")

@app.template_filter('cleanText')
def _jinja2_filter_datetime(text):
    return text.replace("_", " " ).title()


def get_dataset_id_by_name(dataset_name):
    dataset_list = submit(genomics.datasets().list(projectNumber=PROJECT_NUMBER))["datasets"]
    dataset = [x for x in dataset_list if x["name"] == dataset_name]
    if len(dataset) == 1:
        return dataset[0]["id"]
    else:
        return None

def append_metadata(variantSet, metadata_keys = metadata_keys):
        if "metadata" in variantSet:
            for meta in variantSet["metadata"]:
                for k,v in meta.items():
                    if v in metadata_keys:
                        variantSet[meta["key"]] = meta["value"]
        return variantSet

@app.route("/")
def home():
    bq, storage, genomics = connect()
    datasets = submit(genomics.datasets().list(projectNumber=PROJECT_NUMBER))["datasets"]
    return render_template('view.html', **locals())

#==========#
# Datasets #
#==========#

@app.route("/Datasets/")
def Project():
    datasets = submit(genomics.datasets().list(projectNumber=PROJECT_NUMBER))["datasets"]
    section_title = "Datasets"
    breadcrumbs = OrderedDict([("Datasets", "active")])
    return render_template('Datasets.html', **locals())


@app.route('/_add_dataset/', methods = ['POST'])
def add_dataset():
    try:
        if request.form['isPublic'] == "true":
            isPublic = True
        else:
            isPublic = False
        proj_num = submit(genomics.datasets().create(body={"isPublic":isPublic, 
                                           "name":request.form['dataset_name'],
                                           "projectNumber":PROJECT_NUMBER}))
        return jsonify(proj_num)
    except:
        return "ERROR"

@app.route('/_remove_dataset/', methods = ['POST'])
def remove_dataset():
    for i in loads(request.form["datasets"]):
        proj_num = submit(genomics.datasets().delete(datasetId=i))
    return jsonify(proj_num)

#======#
# Jobs #
#======#

@app.route("/Jobs/")
def jobs():
    PROJECT_NUMBER = config["PROJECT_NUMBER"]
    joblist = submit(genomics.jobs().search(body={"projectNumber":PROJECT_NUMBER}))
    joblist = joblist["jobs"]
    breadcrumbs = OrderedDict([("Jobs", "active")])
    return render_template("Jobs.html", **locals())


@app.route("/Datasets/<dataset_name>/")
def Dataset(dataset_name):
    dataset_id = get_dataset_id_by_name(dataset_name)
    if dataset_id is None:
        raise Exception("Multiple databases with the same name not supported yet.")
    variantSets = submit(genomics.variantsets().search(body={"datasetIds":[dataset_id]}))
    variantSets = variantSets["variantSets"]
    # Reorganize variant sets
    variantSets = [append_metadata(x) for x in variantSets]
    breadcrumbs = OrderedDict([("Datasets", url_for('Project')), (dataset_name, "active")])
    return render_template('dataset.html', **locals())

@app.route("/Datasets/<dataset_name>/<variant_set_id>/<ref>/")
def Variants(dataset_name, variant_set_id, ref):
    dataset_id = get_dataset_id_by_name(dataset_name)
    if dataset_id is None:
        raise Exception("Multiple databases with the same name not supported yet.")
    variantSet = submit(genomics.variantsets().search(body={"datasetIds":[dataset_id]}))
    variantSet = variantSet["variantSets"][0]
    variants = submit(genomics.variants().search(body={"referenceName":ref, "variantSetIds":[variant_set_id]}))


    variantSet = append_metadata(variantSet, metadata_keys)

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


    # If a source string is available, use to label the variantSet
    if "source" in variantSet:
        variant_set_name = variantSet["source"]
    else:
        variant_set_name = "variantSet"

    breadcrumbs = OrderedDict([("Datasets", url_for('Project')),
                        (dataset_name, url_for('Dataset', dataset_name=dataset_name)),
                        (variant_set_name,"active"),
                        (ref, "active")])
    return render_template('variant.html', **locals())

@app.route("/<chrom>/<start>/<end>/")
def query(chrom, start, end):
    variants = vcf.Reader(filename="ANNOTATIONS.bcf").fetch(chrom,start,end)
    return render_template('index.html', **locals())


if __name__ == "__main__":
    app.run(debug=True)
    # Enter your Google Developer Project number

