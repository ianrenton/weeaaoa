# "Worked Everything, Everywhere, All At Once" Award Helper Script
# Ian Renton, July 2024
import csv
import json
import math
import os.path
import pathlib
from datetime import datetime

import geopandas as gpd
import pandas as pd
import pyproj
import requests
import shapely
from shapely.geometry import Point, shape

# File to use for raw data caching. If it exists, data will be loaded from it instead of querying the internet. To
# force a redownload, just delete this file.
DATA_FILE = ".cache/data.json"
# File to use for geo data caching. If it exists, data will be loaded from it instead of running the processing. To
# force re-running the processing, just delete this file.
GEO_DATA_FILE = ".cache/geodata.json"
RESULT_FILE = "output.json"

# Static defines
SOTA_ASSOC_BASE_URL = "https://api2.sota.org.uk/api/associations/"
SOTA_REGIONS_BASE_URL = "https://api2.sota.org.uk/api/regions/"
SOTA_ASSOCIATIONS = ["G", "GM", "GW", "GI", "GU", "GJ", "GD"]
BUNKERS_URL = "https://drive.google.com/uc?id=1ea3j9S4VzcDttMPs_9WOj-4L_DzfcUhR"
LIGHTHOUSES_URL = "https://ecaelastats.site/ela_refs.php"
CASTLES_URL = "https://ecaelastats.site/eca_refs.php"
ECL_ECA_FETCH_USER_AGENT = "Mozilla/5.0 (X11; Linux x86_64; rv:123.0) Gecko/20100101 Firefox/123.0"
CIRCLE_TO_POLY_POINTS = 128
SUMMIT_RADIUS_METRES = 250
BUNKER_RADIUS_METRES = 1000
LIGHTHOUSE_RADIUS_METRES = 1000
CASTLE_RADIUS_METRES = 1000

# Static storage
WGS84_TO_OS_GRID_TRANSFORMER = pyproj.Transformer.from_crs(4326, 27700)
OS_GRID_TO_WGS84_TRANSFORMER = pyproj.Transformer.from_crs(27700, 4326)

# Data storage
all_data = []
data_file = pathlib.Path(DATA_FILE)
geo_data_file = pathlib.Path(GEO_DATA_FILE)
result_file = pathlib.Path(RESULT_FILE)

# Main code starts here
# If we don't have an existing data file, start fetching data from the internet.
if not os.path.isfile(DATA_FILE):
    # Retrieve SOTA data
    print("Retrieving SOTA data...")
    for sota_assoc_code in SOTA_ASSOCIATIONS:
        r = requests.get(SOTA_ASSOC_BASE_URL + sota_assoc_code)
        assocs_data = r.json()
        for region in assocs_data["regions"]:
            sota_region_code = region["regionCode"]
            r2 = requests.get(SOTA_REGIONS_BASE_URL + sota_assoc_code + "/" + sota_region_code)
            region_data = r2.json()
            for summit in region_data["summits"]:
                all_data.append({"ref": summit["summitCode"],
                                 "name": summit["name"],
                                 "type": "SUMMIT",
                                 "radiusMetres": SUMMIT_RADIUS_METRES,
                                 "lat": float(summit["latitude"]),
                                 "lon": float(summit["longitude"])})

    # Retrieve castles data
    print("Retrieving ECA data...")
    df = pd.read_html(CASTLES_URL,
                      storage_options={
                          "User-Agent": ECL_ECA_FETCH_USER_AGENT})
    for _, entity in df[0].iterrows():
        all_data.append({"ref": entity.iloc[0],
                         "name": entity.iloc[1],
                         "type": "CASTLE",
                         "radiusMetres": CASTLE_RADIUS_METRES,
                         "lat": float(entity["Latitude"]),
                         "lon": float(entity["Longitude"])})

    # Retrieve lighthouses data
    print("Retrieving ELA data...")
    df = pd.read_html(LIGHTHOUSES_URL,
                      storage_options={
                          "User-Agent": ECL_ECA_FETCH_USER_AGENT})
    for _, entity in df[0].iterrows():
        all_data.append({"ref": entity.iloc[0],
                         "name": entity.iloc[1],
                         "type": "LIGHTHOUSE",
                         "radiusMetres": LIGHTHOUSE_RADIUS_METRES,
                         "lat": float(entity["Latitude"]),
                         "lon": float(entity["Longitude"])})

    # Retrieve bunkers data
    print("Retrieving Bunker data...")
    r = requests.get(BUNKERS_URL, allow_redirects=True)
    reader = csv.DictReader(r.text.split('\n'))
    for entity in reader:
        all_data.append({"ref": entity["UKBOTA Reference"],
                         "name": entity["Bunker Name"],
                         "type": "BUNKER",
                         "radiusMetres": BUNKER_RADIUS_METRES,
                         "lat": float(entity["Latitude"]),
                         "lon": float(entity["Longitude"])})

    # Write cache file
    print("Writing cache file...")
    data_file.parent.mkdir(exist_ok=True, parents=True)
    data_file.write_text(json.dumps(all_data, indent=2))

else:
    print("Reading cache file...")
    with open(data_file) as f:
        all_data = json.load(f)

print(str(len(all_data)) + " entities found.")

print("Converting data for GeoPandas...")
start = datetime.now()

# Convert all lat/lons to OS grid reference. This will break totally for data outside the UK, but it saves a lot of
# hassle in GeoPandas because everything can be in metres.
for entity in all_data:
    os_grid_ref = WGS84_TO_OS_GRID_TRANSFORMER.transform(entity["lat"], entity["lon"])
    entity["northing"] = os_grid_ref[0]
    entity["easting"] = os_grid_ref[1]

# Prepare a polygons for each entity in OS grid reference space. We will need this later whether or not we have a
# cached geo data file containing the overlap segments.
for entity in all_data:
    gs = gpd.GeoSeries(Point(entity["northing"], entity["easting"]), crs=27700)
    buffer_gs = gs.buffer(entity["radiusMetres"], CIRCLE_TO_POLY_POINTS)
    entity["polygon"] = buffer_gs[0]

all_buffers_geoseries = gpd.GeoSeries(list(map(lambda p: p["polygon"], all_data)))

# Now run the processing, if we don't already have the result data
if not os.path.isfile(geo_data_file):
    # Assemble a GeoDataFrame containing all the entities.
    data = {'name': list(map(lambda p: p["ref"] + " " + p["name"], all_data)),
            'id': range(0, len(all_data)),
            'geom': all_buffers_geoseries}

    df = pd.DataFrame(data, columns=['name', 'id', 'geom'])
    gdf = gpd.GeoDataFrame(df, geometry='geom', crs=27700)

    print("Finding overlaps...")
    # Code from https://gis.stackexchange.com/questions/387773/count-overlapping-features-using-geopandas
    buffer_size = 0.1
    bounds = gdf.geometry.convex_hull.exterior.buffer(buffer_size).unary_union
    new_polys = list(shapely.ops.polygonize(bounds))
    # Removing the full merged polygons (first is always index 0,
    # subsequent will be the first of their own 'bunches' identified as disjoint from other 'bunches')
    bad_poly_idx = [0]
    while new_polys[max(bad_poly_idx)].disjoint(new_polys[-1]):
        for idx in range(max(bad_poly_idx), len(new_polys)):
            if new_polys[max(bad_poly_idx)].disjoint(new_polys[idx]):
                bad_poly_idx += [idx]
                break
    new_polys = [new_polys[i].buffer(-buffer_size) for i in range(len(new_polys)) if i not in bad_poly_idx]
    # count layers and track IDs of overlapping features
    gdf_with_overlap_polys = gpd.GeoDataFrame(geometry=new_polys)
    gdf_with_overlap_polys['layers'] = sum(
        [gdf_with_overlap_polys.geometry.intersects(poly) for poly in gdf.geometry.buffer(buffer_size).values])
    gdf_with_overlap_polys['piece'] = gdf_with_overlap_polys.index

    runtime = datetime.now() - start
    print("Generated " + str(len(gdf_with_overlap_polys.index)) + " overlap polys in " + str(
        runtime.total_seconds()) + " seconds.")

    print("Writing cache file...")
    geo_data_file.parent.mkdir(exist_ok=True, parents=True)
    geo_data_file.write_text(gdf_with_overlap_polys.to_json())
    print("Written " + str(int(os.path.getsize(GEO_DATA_FILE) / 1024 / 1024)) + " MB.")

else:
    print("Reading cache file...")
    with open(geo_data_file) as f:
        gdf_with_overlap_polys = gpd.read_file(f)

# Now iterate over the overlap poly features. Only care about ones with more than one "layer" (i.e. overlapping original
# entity). For each such overlap poly, create a test point inside it, and see which original entities it's in range of.
# Store that example point along with the list of entities in range.
print("Getting entity lists for overlap polygons...")
start = datetime.now()
overlap_data = []
for feature in gdf_with_overlap_polys.iterfeatures():
    if feature["properties"]["layers"] > 1:
        test_point = shape(feature["geometry"]).representative_point()
        overlapping_entity_names = []
        for test_entity in all_data:
            north_dist = abs(test_point.x - test_entity["northing"])
            east_dist = abs(test_point.y - test_entity["easting"])
            dist = math.sqrt(east_dist * east_dist + north_dist * north_dist)
            if dist < test_entity["radiusMetres"]:
                overlapping_entity_names.append(test_entity["ref"] + " " + test_entity["name"])
        test_point_wgs84 = OS_GRID_TO_WGS84_TRANSFORMER.transform(test_point.x, test_point.y)
        overlap_data.append({"lat": test_point_wgs84[0],
                             "lon": test_point_wgs84[1],
                             "entities": overlapping_entity_names})

runtime = datetime.now() - start
print("Assessed " + str(len(gdf_with_overlap_polys.index)) + " overlap polys in " + str(
    runtime.total_seconds()) + " seconds.\n")

# Sort overlap data by number of entities
overlap_data_sorted = sorted(overlap_data, key=lambda p: len(p["entities"]), reverse=True)
print("Top 5 regions by overlapping entity count:")
for i in range(0, 4):
    print(str(overlap_data_sorted[i]["lat"]) + " " + str(overlap_data_sorted[i]["lon"]) + " " + str(
        len(overlap_data_sorted[i]["entities"])) + " " + str(overlap_data_sorted[i]["entities"]))

print("\nWriting results file...")
result_file.parent.mkdir(exist_ok=True, parents=True)
result_file.write_text(json.dumps(overlap_data_sorted, indent=2))
print("Done.")
