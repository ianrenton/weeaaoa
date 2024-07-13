# Provide the index of a particular region in the output data, get a list formatted for the Keene State College Polyline
# Tool

import json
import pathlib

RESULT_FILE = "output.json"
INDEX = 0

result_file = pathlib.Path(RESULT_FILE)
with open(result_file) as f:
    all_data = json.load(f)
    poly = all_data[INDEX]["polygon"]
    for latlon in poly:
        print(str(latlon[1]) + "," + str(latlon[0]))