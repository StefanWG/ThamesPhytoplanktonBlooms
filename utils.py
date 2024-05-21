'''
Various helper functions for smoothing and data processing
'''
import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
from copy import deepcopy
from sympy import sympify, symbols
from PIL import Image

#TODO: Create function for weeksSinceStat and dates

##############################
# Data preprocessing
##############################

def readData(loc):
    data = pd.read_csv(f"data/processed/{loc}.csv")
    data["date"] = pd.to_datetime(data["date"])
    return data

def processData(filePath, vars = ["date","temp","flow","SRP","sun3d","sun","silicon","nitrate",
                "chl", "diatoms", "p-chloro","n-chloro","totCyano", "totCrypto", "totPhyto"],
                outputFolder = "data/processed"):
    
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
        data.drop(columns=["index"], inplace=True)

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
        rm = readData("runnymede")
        weeksSinceStart = [((i - rm["date"][0]).days + 1) // 7 + 1 for i in rm["date"]]
        rm["weeks"] = weeksSinceStart

        data["weeks"] = [((i - rm["date"][0]).days + 1) // 7 + 1 for i in data["date"]]
        rm = rm[["weeks", "sun", "sun3d"]]

        data = data.merge(rm, on="weeks", how="left")

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


def createCSVDoc(outputFile = "data/processed/documentation.csv"):
    '''
    Create documentation for CSVs.

    Inputs:
        - outputFile (str): file path for output 
    '''
    d = [
        {"name":"date", "fullName": "Date"},
        {"name":"temp", "fullName": "Temperature", "unit":"°C"},
        {"name":"flow", "fullName": "River Flow", "unit":'m3/s'},
        {"name":"SRP", "fullName": "Soluble Reactive Phosphorus", "unit":"µg/l"},
        {"name":"sun3d", "fullName": "Sunlight (3 Day Mean)", "unit":"Hours"},
        {"name":"sun", "fullName": "Sunlight", "unit":"Hours"},
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

##############################
# Smoothing
##############################

def convolve1d(vec, filter):
    '''
    Creates a convoluted vector given a specific filter.
    Fills vector with 0s on both ends so len(vec) == len(vOut).

    Inputs:
        - vec (array): vector to convolve
        - fitler (array): convolution filter

    Outputs:
        - vOut (array) convolved vector
    '''
    mid = len(filter) // 2
    # Create output array
    vOut = [0 for i in range(len(vec))]

    for i in range(len(vec)):
        # Sub vector for convolution
        vFilt = np.array(vec[(max(0, i - mid)):(min(len(vec)-1, i+mid+1))])
        # Pad vector with 0s if too short
        if len(vFilt) < len(filter):
            if i < mid:
                vFilt = np.append(np.zeros(len(filter) - len(vFilt)), vFilt)
            else:
                vFilt = np.append(vFilt, np.zeros(len(filter) - len(vFilt)))
        # Perform convolution
        mask = vFilt != 0
        vOut[i] = sum(vFilt * filter) / sum(mask * filter)
    return vOut

def smooth(vec,sigma=None, type="gaussian"):
    '''
    Smooth and gap fill time series data with desired filter.
    Gaps must be 0.0 and not np.nan.

    Inputs:
        - vec (array): time series data
        - sigma (int): sigma for gaussian curve (default: None)
        - type (str): type of filter to use for convolution; either
                      "gaussian" or "uniform" (default: "gaussian)

    Outputs:
        - smoothed (array): smoothed and gapfilled time series
    '''
    # Create Gaussian Curve
    x = np.arange(-3*sigma,3*sigma+1)
    if type == "gaussian":
        assert sigma is not None
        filter = np.exp((-(x/sigma)**2)/2.0)
    elif type == "uniform":
        filter = [1.0 for i in x]
    else:
        raise ValueError("Invalid filter type")
    # Perform Smoothing and Gapfilling
    smoothed = convolve1d(vec, filter)
    return smoothed

def getSmoothedData(dates, vec, sigma=None, type="gaussian"):
    '''
    Get smoothed and gap filled data for column in data table

    Inputs:
        - dates (array): list of dates corresponding to time series data
        - vec (array): time series data
        - sigma (int): sigma for gaussian curve (default: None)
        - type (str): type of filter to use for convolution; either
                      "gaussian" or "uniform" (default: "gaussian)

    Outputs:
        - smoothed (array): smoothed and gapfilled data
    '''
    # Create array with entry for each date, dates with gaps have 0
    # day0 = dates.iloc[:1].values[0]
    weeksSinceStart = [((i - dates[0]).days + 1) // 7 + 1 for i in dates]
    vMap = {a:b for a,b, in zip(weeksSinceStart, vec)}
    toSmooth = np.array([vMap[i+1] if i+1 in vMap else 0 for i in range(max(weeksSinceStart))], dtype=float)
    toSmooth[np.isnan(toSmooth)] = 0
    # Smooth
    smoothed = np.array(smooth(toSmooth, sigma, type))
    return smoothed


##############################
# Finding blooms and cessations
##############################

def getThresholds(vec, peakMins, peakMaxs, weeksSinceStart, confInterval=90):
    '''
    Get thresholds during which blooms occur 

    Inputs:
        - vec (array): time series data to get thresholds from
        - peakMins (array): array of starting point of peaks
        - peakMaxes (array): array of ending point of peaks
        - weeksSinceStart (array): array of weeks passed since first data point
                                   for all data points
        - confInterval (float): only increase confInterval % of values when 
                                computing threshold (default: 90)

    Outputs:
        - thresholds (tuple): minimum threshold and maximum threshold (both floats)
    '''
    # Get indices of peaks
    peakIdxs = np.concatenate([np.arange(i1, i2+1) for i1, i2, in zip(peakMins, peakMaxs)])
    # Create array to extract thresholds from
    varMap = {a:b for a,b, in zip(weeksSinceStart, vec)}    
    var = np.array([varMap[i+1] if i+1 in varMap else np.nan for i in range(max(weeksSinceStart))])
    thresholds = (np.nanpercentile(var[peakIdxs], (100-confInterval) / 2), 
                  np.nanpercentile(var[peakIdxs], (100-confInterval) / 2 + confInterval))
    return thresholds

def getBlooms(dates, vec, sigma=None, type="gaussian"):
    '''
    Get dates during with blooms and cessations occured.
    This is given by peaks using the scipy find_peaks algorithm.

    Inputs:
        - dates (array): list of dates corresponding to time series data
        - vec (array): time series data
        - sigma (int): sigma for gaussian curve (default: None)
        - type (str): type of filter to use for convolution; either
                      "gaussian" or "uniform" (default: "gaussian)

    Outputs:
        - smoothed (array): array of smoothed and gap filled time series data
        - peaks (array): array of peaks (this corresponds to where bloom becomes cessation)
        - starts (array): array of bloom starting points
        - ends (array): array of cessation ending points


    '''
    smoothed = getSmoothedData(dates, vec, sigma, type)
    peaks, _ = scipy.signal.find_peaks(smoothed, prominence=np.max(smoothed) /20) # TODO: Modify this based on time series data
    widths = scipy.signal.peak_widths(smoothed, peaks)
    starts = np.round(widths[2]).astype(int)
    ends = np.round(widths[3]).astype(int)

    return smoothed, peaks, starts, ends

##############################
# Utility functions
##############################

def rmse(pred, obs):
    '''
    Calculate root mean squared error (RMSE)
    Note: both must be numpy arrays

    Inputs:
        - pred (np.array): array of predicted values
        - obs (np.array): array of observed values
    
    Outputs:
        - rmse (float): rmse between predicted and observed
    '''
    assert len(pred) == len(obs)
    rmse = np.sqrt(np.mean((pred - obs)**2))
    return rmse

def getInflectionPoints(vec):
    '''
    Get inflection points using the second derivative.
    Since points are at discrete x-values (as opposed to 
    a continous function), inflection points 
    are where the line crosses the x-axis.

    Inputs:
        - vec (array): array from which the inflection points 
                    are extracted
    
    Output:
        - infls (array): array of indices of inflection points
    '''
    secondDeriv = np.diff(vec, 2)
    infls = np.where(np.diff(np.sign(secondDeriv)))[0] 

    return infls

def groupPlots():
    '''
    Create image for each type of plot with subplots for each station
    '''
    # TODO: Not robust and poor comment

    for var in ["chl", "diatoms", "pchloro", "nchloro", "totCyano", "totPhyto", "linModel"]:
        fig = plt.figure(figsize=(20,20))

        for i, loc in enumerate(["cherwell", "kennet", "ock", "ray", "runnymede","thame"]):
            ax = fig.add_subplot(2,3,i+1)
            img = np.asarray(Image.open(f'figures/{loc}_{var}.jpg'))
            ax.axis('off')
            ax.imshow(img)
            ax.set_title(loc)

        fig.tight_layout()
        fig.subplots_adjust(bottom=0, top=1, hspace= -0.68 if var != "linModel" else -.5)
        fig.savefig(f"figures/{var}.jpg", bbox_inches="tight")
