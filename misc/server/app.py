import os
import json
import io
from flask import Flask, render_template, request, redirect, jsonify, stream_with_context, Response
from hobotnica import *
import sys
from base64 import b64encode

app = Flask(__name__)

@app.route('/predict', methods=['GET', 'POST'])
def predict():
    if request.method == 'POST':
        countmatrix = request.files['countmatrix']
        annotation = request.files['annotation']

        print("=== SIZE ===")
        t1 = len(countmatrix.read())
        t2 = len(annotation.read())
        print(t1, t2)
        print("============")
        countmatrix.seek(0)
        annotation.seek(0)

        hbt_container = create_hobotnica(countmatrix, annotation, True)
        try:
            outstream = call_hobotnica(
                hbt_container,
                debug=True
            )
        except:
            print("\n=== FAILED TO PROCESS HOBOTNICA ===\n")
            outstream = "Failed"

        return Response(outstream)
    else:
        return 'No files specified!'

