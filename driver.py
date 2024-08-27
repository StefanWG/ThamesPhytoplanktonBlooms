import argparse 
from dataProcessing import processData, createCSVDoc

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--csvs", dest = "csvs", default=False, help="Process CSVs", action="store_true")

args = parser.parse_args()


if args.csvs:
    for location in ["runnymede", "ock", "kennet", "thame", "cherwell", "ray"]:
        print(f"Processing data for {location}...")
        processData(f"data/{location}.csv")

    createCSVDoc()

