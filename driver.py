import subprocess 
from utils import *
import argparse 
from linModel import LinModel

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--csvs", dest = "csvs", default=False, help="Process CSVs", action="store_true")
parser.add_argument("-p", "--plots", dest = "plots", default=False, help="Generate Plots", action="store_true")
parser.add_argument("-m", "--models", dest = "models", default=False, help="Run Models", action="store_true")

args = parser.parse_args()


if args.csvs:
    for location in ["ray", "ock", "kennet", "thame", "cherwell", "runnymede"]:
        print(f"Processing data for {location}...")
        processData(f"data/{location}.csv")

    createCSVDoc()

if args.plots:
    for location in ["ray", "ock", "kennet", "thame", "cherwell", "runnymede"]:
        print(f"Making plots for {location}...")
        subprocess.run(["python3", "plots.py", f"figures/{location}_", location])
 
if args.models:
    for location in ["ray", "ock", "kennet", "thame", "cherwell", "runnymede"]:
        print(f"Running for {location}...")
        model = LinModel(location)
        model.plot(output=f"figures/{location}_linModel.jpg")
