import pandas as pd
import numpy as np
from utils import *
from smoothing import *
import os 
import datetime 

def readData(loc):
    '''
    Read data from a CSV file.

    Inputs:
        - loc (str): Location of the CSV file.

    Returns:
        - data (pandas.DataFrame): The loaded data.
    '''
    data = pd.read_csv(f"data/processed/{loc}.csv")
    data["date"] = pd.to_datetime(data["date"])
    return data

def processData(loc, filePath, vars = COLS,
                outputFolder = "data/processed"):
    '''
    Process data from a CSV file.

    Inputs:
        - loc (str): Location of the CSV file.
        - filePath (str): Path to the CSV file.
        - vars (list): List of variables to include in the processed data. Default is COLS.
        - outputFolder (str): Output folder for the processed data. Default is "data/processed".

    Returns:
        None
    '''
    assert outputFolder != "data" # This will override files
    if outputFolder[-1] == "/":
        outputFolder = outputFolder[:-1]

    data = pd.read_csv(filePath, low_memory=False)

    if filePath in ["data/thamesData.csv", "data/runnymede.csv"]:
        data.rename(columns = {"Date":"date",
                        "Water temperature":"temp",
                        "Chlorophyll-a (_g/L)":"chl",
                        "Mean daily flow (m3/s)":"flow",
                        "Sunshine duration (3-d average)":"sun3d",
                        "Sunshine duration ":"sun",
                        "Total phytoplankton":"totPhyto",
                        "Diatoms":"diatoms",
                        "Nano-chlorophytes":"n-chloro",
                        "Pico-chlorophytes":"p-chloro",
                        "Total Cryptophytes":"totCrypto",
                        "Total Cyanobacteria":"totCyano",
                        "Soluble Reactive P":"SRP",
                        "Silicon (mg-Si/L)":"silicon",
                        "Dissolved nitrate (NO3)":"nitrate"
                        }, inplace=True)
        # Fix dates
        m = [str(i).split(".")[0] if not np.isnan(i) else "-" for i in data["Month"]]
        d = [str(i).split(".")[0] if not np.isnan(i) else "-" for i in data["Day"]]
        y = [str(i).split(".")[0][-2:] if not np.isnan(i) else "-" for i in data["Year"]]
        data["date"] = [pd.to_datetime(f"{a}/{b}/20{c}") if not any([i == "-" for i in [a,b,c]]) else np.nan for a, b, c, in zip(m, d, y)]
        data = data[data["temp"] > -900]
        data.reset_index(inplace=True)
        data.drop(columns=["index", "sun", "sun3d"], inplace=True)

         # Sun
        closestStations = getClosestSunStations()
        sunStation = closestStations[closestStations["label"] == loc].iloc[0]["closestStation"]
        sunData = getSunData(sunStation)
        sunData["date"] = pd.to_datetime(sunData["date"])

        data = pd.merge(data, sunData, on="date", how="left")
        data = data.reset_index().drop(columns=["index"])
        data_obj = data.select_dtypes('object')

    else:
        # Clean dates
        data = data[~pd.isna(data["Day"])]
        m = [str(i).split(".")[0] if not np.isnan(i) else "-" for i in data["Month"]]
        d = [str(i).split(".")[0] if not np.isnan(i) else "-" for i in data["Day"]]
        y = [str(i).split(".")[0][-2:] if not np.isnan(i) else "-" for i in data["Year"]]
        y = ["0" + i if len(i) == 1 else i for i in y]
        data["date"] = [pd.to_datetime(f"{a}/{b}/20{c}") if not any([i == "-" for i in [a,b,c]]) else np.nan for a, b, c, in zip(m, d, y)]

        data["Temperature (oC)"] = [np.nan if d == " " or pd.isna(d) else d for d in data["Temperature (oC)"]]

        data.rename(columns={
            "Mean daily flow (m3/s)":"flow",
            "Temperature (oC)":"temp",
            'Chlorophyll-a (μg/L)':"chl",
            "Soluble reactive phosphorus (μg/L)":"SRP",
            "Dissolved silicon (mg Si/L) ":"silicon",
            "Dissolved nitrate (NO3)":"nitrate",
            "Green_diatoms":"diatoms",
            "Green_pico":"p-chloro"
        }, inplace=True)

        # Calculate cell counts
        data["n-chloro"] = [a+b for a,b in zip(data["Green_meso"], data["Green meso 2"])]
        data["totCrypto"] = [a+b+c for a,b,c in zip(data["Crypto 1"], data["Crypto 2"], data["Crypto 3"])]
        data["totCyano"] = [a+b+c+d for a,b,c,d, in zip(data["Cyano1"], data["Cyano2"], data["Cyano3"], data["Cyano4"])]
        data["totPhyto"] = [a+b+c+d+e for a,b,c,d,e in zip(data["diatoms"], data["totCyano"], data["n-chloro"], data["p-chloro"], data["totCrypto"])]

        # Sun
        closestStations = getClosestSunStations()
        sunStation = closestStations[closestStations["label"] == loc].iloc[0]["closestStation"]
        sunData = getSunData(sunStation)
        sunData["date"] = pd.to_datetime(sunData["date"])

        data = pd.merge(data, sunData, on="date", how="left")
        data = data.reset_index().drop(columns=["index"])
        data_obj = data.select_dtypes('object')

        data[data_obj.columns] = data_obj.apply(lambda x: x.str.strip())
        data = data.replace({"":np.nan})

        # Cut off data when a single time series ends
        s,e  = 0, len(data)
        for var in vars:
            if var in ["date", "sun", "sun3d"]:
                continue
            exists = np.where(~np.isnan(np.array(data[var].astype(float))))
            s = max(s, exists[0][0])
            e = min(e, exists[0][-1])
        data = data.iloc[s:e+1]

    outputFp = f'{outputFolder}/{filePath.split("/")[-1]}'

    data = data[vars]


    data.to_csv(outputFp, index=False)
    print(f"CSV saved to {outputFp}.")

def getSunData(station):
    '''
    Get solar radiation data for a specific station.

    Inputs:
        - station (str): Station name.

    Returns:
        - dailySR (pandas.DataFrame): Solar radiation data.
    '''
    files = os.listdir(f"data/sun/{station}")
    files = [f for f in files if f != ".DS_Store"]

    df = pd.read_csv(f"data/sun/{station}/{files[0]}", skiprows=75)
    df = df[df["ob_hour_count"] == 1]
    df["date"] = pd.to_datetime(df["ob_end_time"])
    df["date"]  =df["date"].apply(lambda x : datetime.datetime.strftime(x, "%Y-%m-%d"))
    dailySR = df.groupby(by="date", as_index=False).agg({"glbl_irad_amt":"mean"})

    for f in files[1:]:
        df = pd.read_csv(f"data/sun/{station}/{f}", skiprows=75)
        df = df[df["ob_hour_count"] == 1]
        df["date"] = pd.to_datetime(df["ob_end_time"])
        df["date"]  =df["date"].apply(lambda x : datetime.datetime.strftime(x, "%Y-%m-%d"))
        d = df.groupby(by="date", as_index=False).agg({"glbl_irad_amt":"mean"})
        dailySR = pd.concat([dailySR, d])


    dailySR = dailySR.rename(columns={"glbl_irad_amt":"sun"})
    dailySR = dailySR.sort_values(by="date")
    dailySR["sun3d"] = dailySR["sun"].rolling(3).mean()
    return dailySR

def smoothDataset(filePath):
    '''
    Smooth the dataset.

    Inputs:
        - filePath (str): Path to the dataset CSV file.

    Returns:
        None
    '''
    data = pd.read_csv(filePath)
    data["date"] = pd.to_datetime(data["date"])

    d = {k:getSmoothedData(data['date'], data[k], SIGMA) for k in data.columns[1:]}
    df = pd.DataFrame(data=d)
    df["date"] = getGapFilledDates(getWeeksSinceStart(data["date"]), data["date"][0])

    df.to_csv(filePath, index=False)


def createCSVDoc(outputFile = "data/processed/documentation.csv"):
    '''
    Create documentation for CSVs.

    Inputs:
        - outputFile (str): File path for the output CSV.

    Returns:
        None
    '''
    d = [
        {"name":"date", "fullName": "Date"},
        {"name":"temp", "fullName": "Temperature", "unit":"°C"},
        {"name":"flow", "fullName": "River Flow", "unit":'m3/s'},
        {"name":"SRP", "fullName": "Soluble Reactive Phosphorus", "unit":"µg/l"},
        {"name":"sun3d", "fullName": "Solar Radiation (3 Day Mean)", "unit":"kJ/h"},
        {"name":"sun", "fullName": "Solar Radition", "unit":"kJ/h"},
        {"name":"silicon", "fullName": "Dissolved Silicon (SI)", "unit":"mg/l"},
        {"name":"nitrate", "fullName": "Dissolved Nitrate (NO3)", "unit":"mg/l"},
        {"name":"chl", "fullName": "Chlorophyll-a", "unit":"µg/l"},
        {"name":"diatoms", "fullName": "Diatoms", "unit":"Cell Count"},
        {"name":"p-chloro", "fullName": "Pico-chlorophytes", "unit":"Cell Count"},
        {"name":"n-chloro", "fullName": "Nano-chlorophytes", "unit":"Cell Count"},
        {"name":"totCyano", "fullName": "Total Cyanobacteria", "unit":"Cell Count"},
        {"name":"totCrypto", "fullName": "Total Cryptophytes", "unit":"Cell Count"},
        {"name":"totPhyto", "fullName": "Total Phytoplankton", "unit":"Cell Count"},
    ]

    pd.DataFrame.from_dict(d).to_csv(outputFile, index=False)

    