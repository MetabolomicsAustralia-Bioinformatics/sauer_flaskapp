from flask import Flask, render_template, url_for, request, redirect, make_response, flash
from flask_wtf import FlaskForm
#from flask_sqlalchemy import SQLAlchemy
from datetime import datetime

import io
import csv
import numpy as np
import pandas as pd

import SauerFunction as sf
from SauerClass import Record, Labelling

app = Flask(__name__)
app.config["SECRET_KEY"] = "xdVag8zTLPMv4UFMyanv"

# Set globals
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///test.db'
ALLOWED_EXTENSIONS = {'csv'}
#db = SQLAlchemy(app)


def make_response_string(labelling_list, max_results_length):
    """Make a string to be passed to make_response(), eventually to be written out as .csv.
    
    PARAMS
    ------
    labelling_list: list of label class objects; one per species.
    
    RETURNS
    -------
    stream_out: str; comma-separated and \n newline-delimited string to be written out.
    """
    # Make header row
    stream_out = "Species,Labelling Source,Sample Name,Labelling %,"
    for i in range(max_results_length-1):
        stream_out += "m" + str(i) + ","
    stream_out += "\n"

    # Make contents
    for label in labelling_list:
        species = label.get_species()
        label_dict = label.get_label_dict()
        names = label_dict.keys()
        names = sorted(names)
        i = 0
        for name in names:
            for key, value in label_dict.items():
                if name == key:
                    stream_out += species + ','
                    if len(key.split(',')) == 2:
                        stream_out += key.split(',')[1].strip('"') + ','
                        stream_out += key.split(',')[0].strip('"') + ','
                    else:
                        stream_out += ' ,'
                        stream_out += key + ','
                    for val in value:
                        stream_out += str(val) + ', '
                    stream_out += '\n'
    
    # Remove all equals-signs
    stream_out = stream_out.replace("=", ",")
    
    return stream_out


@app.route("/", methods=['POST', 'GET'])
def index():
	return render_template("index.html")


@app.route("/download", methods=["GET"])
def download():
    print("Download button clicked.")
    return redirect("/")


@app.route('/transform', methods=["GET", "POST"])
def transform_view():
    f = request.files['input_file1']
    if not f:
        return "Main input file not found!"
    g = request.files['input_file2']
    if not f:
        return "Formulae input file not found!"

    stream = io.StringIO(f.stream.read().decode("UTF8"), newline=None)
    csv_input = csv.reader(stream)
    contents = []
    for row in csv_input:
        contents.append(row)

    record_ls = sf.read_hunter_contents(contents)
    print("main input file uploaded!")

    stream = io.StringIO(g.stream.read().decode("UTF8"), newline=None)
    csv_input = csv.reader(stream)
    contents2 = []
    for row in csv_input:
        contents2.append(row)

    atomic_composition, N_dict = sf.read_atomic_contents(contents2)
    print("Formulae file uploaded!")

    labelling_list = []
    max_results_length = 0
    for record in record_ls:
        results_dict = sf.calculate_labelling(record, N_dict, atomic_composition)
        for key, value in results_dict.items():
            if len(value) > max_results_length:
                max_results_length = len(value)
        labelling_list.append(Labelling(record.get_name(), results_dict))

    stream_out = make_response_string(labelling_list, max_results_length)
    print(stream_out)

    response = make_response(stream_out)
    response.headers["Content-Disposition"] = "attachment; filename=result.csv"

    return response


if __name__ == "__main__":
	app.run(debug=True)
