import numpy as np
import pandas as pd

LOCATIONS = ["runnymede", "ock", "cherwell", "ray", "kennet", "thame"]
TARGETVARS = ["diatoms", "p-chloro", "n-chloro", "totCyano", "totCrypto", "totPhyto", "chl"]
INDVARS = ["temp", "flow", "sun", "sun3d", "SRP", "silicon", "nitrate"]
COLS = ["date"] + INDVARS + TARGETVARS
BLOOMMETHODS = ["find_peaks", "kneedle", "inflection"]
SIGMA = 3

LABELS = {
    "diatoms":"Diatoms", 
    "n-chloro":"Nano-chlorophytes", 
    "p-chloro":"Pico-chlorophytes", 
    "totCyano":"Cyanobacteria", 
    "totCrypto":"Total Cryptophytes",
    "totPhyto":"Total Phytoplankton",
    "chl":"Chlorophyll"
}

LABELSUNITS = {
    "diatoms":"Diatom Cells\n(m/l)", 
    "temp":"Water Temperature\n(°C)",
    "flow":r'River Flow' "\n" r'(m$^3$/s)',
    "sun3d":'Sunlight Hours\n3 Day Mean)',
    "n-chloro":"Nano-chlorophytes\n(m/l)", 
    "p-chloro":"Pico-chlorophytes\n(m/l)",
    "totCyano":"Cyanobacteria\n(m/l)",
    "totPhyto":"Total Phytoplankton\n(m/l)",
    "totCrypto":"Total Cryptophytes\n(m/l)",
    "chl":"Chlorophyll"
}

# Color Palattes
COLORS = ["b","r", "g", "y", "k", "m", "c"]
COLORSTAB = ["tab:brown", "black", "tab:blue", "tab:olive"]

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

def getWeeksSinceStart(dates):
    '''
    Get number of weeks since first data point for all
    data points (indexing begins at 1).

    Inputs:
        - dates (array): array of dates
    
    Outputs:
        - w (array): list of weeks since start
    '''
    w = [((i - dates[0]).days + 1) // 7 + 1 for i in dates]
    return w

def getGapFilledDates(weeksSinceStart, startDate):
    '''
    Get list of dates that has been gap filled.

    Inputs:
        - weeksSinceStart (array): list of weeks since start
                                   for each data point
        - startDate (date): date of first data point 

    Outputs:
        - dates (array): gap filled array of dates
    '''
    dates = [startDate + pd.Timedelta(weeks=i) for i in range(max(weeksSinceStart))]
    return dates
