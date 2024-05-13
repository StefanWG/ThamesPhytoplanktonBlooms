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

#TODO: Create function for weeksSinceStat and dates

##############################
# Data preprocessing
##############################

def readData(file="data/thamesData.csv"):
    '''
    Read in data from csv, change column names,
    and process dates

    Inputs:
        - file (str): filepath to csv (default: "data/thamesData.csv")

    Outputs:
        - date (pd.Dataframe): pandas dataframe of processed csv 
    '''
    data = pd.read_csv(file)

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

    return data

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
    weeksSinceStart = [((i - dates[0]).days + 1) // 7 + 1 for i in dates]
    vMap = {a:b for a,b, in zip(weeksSinceStart, vec)}
    toSmooth = np.array([vMap[i+1] if i+1 in vMap else 0 for i in range(max(weeksSinceStart))])
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
    peaks, _ = scipy.signal.find_peaks(smoothed, prominence=1200) # TODO: Modify this based on time series data
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
    # TODO: Check for numpy, and implement for list
    rmse = np.sqrt(np.mean((pred - obs)**2))
    return rmse