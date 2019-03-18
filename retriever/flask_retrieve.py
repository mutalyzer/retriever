from flask import Flask, jsonify
from flask import request

from retriever.parsers.gff_other import get_raw_record

app = Flask(__name__)


@app.route('/get-feature')
def hello_world():
    feature_id = request.args.get('id', default='NM_002001.3', type=str)
    model = get_raw_record(feature_id)
    return jsonify(results=model, sort_keys=True)


@app.route('/get-raw-model')
def my_route():
    feature_id = request.args.get('id', default='NM_002001.3', type=str)
    model = get_raw_record(feature_id)
    return jsonify(results=model)
