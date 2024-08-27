import numpy as np
import pandas as pd
import haversine as hs
import datetime
from enum import Enum

LOCATIONS = ["runnymede", "ock", "cherwell", "ray", "kennet", "thame"]
TARGETVARS = ["diatoms", "p-chloro", "n-chloro", "totCyano","chl"]
INDVARS = ["temp", "flow", "sun", "sun3d", "SRP", "silicon", "nitrate"]
COLS = ["date"] + INDVARS + TARGETVARS
BLOOMMETHODS = ["find_peaks", "kneedle"]
SIGMA = 3

LOCATIONMAP = {
    "runnymede":"River Thames at Runnymede",
    "ock":"River Ock",
    "ray":"River Ray",
    "kennet":"River Kennet",
    "thame":"River Thame",
    "cherwell":"River Cherwell",
    "na":""
}

LABELS = {
    "diatoms":"Diatoms", 
    "n-chloro":"Nano-chlorophytes", 
    "p-chloro":"Pico-chlorophytes", 
    "totCyano":"Cyanobacteria", 
    "totCrypto":"Total Cryptophytes",
    "totPhyto":"Total Phytoplankton",
    "chl":"Chlorophyll",
    "flow":"River Flow",
    "sun3d":"Solar Radiation \n(3-Day Mean)",
    "temp":"Temp",
    "SRP":"Phosphorus",
    "silicon":"Silicon",
    "nitrate":"Nitrate",
    "t_q":"TQ"
}

LABELSUNITS = {
    "diatoms":"Diatom Cells\n(cells / mL)", 
    "temp":"Water Temperature\n(°C)",
    "flow":r'River Flow' "\n" r'(m$^3$/s)',
    "sun3d":'Solar Radiation\n(3 Day Mean, kJ/h)',
    "n-chloro":"Nano-chlorophytes\n(cells / mL)", 
    "p-chloro":"Pico-chlorophytes\n(cells / mL)",
    "totCyano":"Cyanobacteria\n(cells / mL)",
    "totPhyto":"Total Phytoplankton\n(cells / mL)",
    "totCrypto":"Total Cryptophytes\n(cells / mL)",
    "chl":"Chlorophyll",
    "SRP":"Soluble Reactive Phosphorus\n(µg/L)",
    "silicon":"Silicon\n(mg/L)",
    "nitrate":"Nitrate\n(mg/L)"
}

# Color Palattes
COLORS = ["b","r", "g", "y", "k", "m", "c"]
COLORSTAB = ["tab:brown", "black", "tab:blue", "tab:olive"]
COLORSXKCD = ["xkcd:snot green", "xkcd:dusty blue", "xkcd:pumpkin"]

def rmse(pred: np.ndarray, obs: np.ndarray) -> float:
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

def getWeeksSinceStart(dates: np.ndarray) -> list:
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

def getGapFilledDates(weeksSinceStart:np.ndarray, startDate: datetime) -> list:
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

def fStat(obs: np.ndarray, pred: np.ndarray) -> float:
    '''
    Calculate F-statistic for two binary arrays.

    Inputs:
        - obs (np.array): array of observed values (binary)
        - pred (np.array): array of predicted values (binary)
    
    Outputs:
        - f (float): F-statistic
    '''
    overlap = obs & pred
    f = sum(overlap) / (sum(obs) + sum(pred) - sum(overlap))
    return f

def NSE(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    '''
    Calculate Nash-Sutcliffe Efficiency (NSE) for two arrays.

    Inputs:
        - y_true (np.array): array of true values
        - y_pred (np.array): array of predicted values
    
    Outputs:
        - nse (float): Nash-Sutcliffe Efficiency
    '''
    return 1 - np.sum((y_true - y_pred)**2) / np.sum((y_true - np.mean(y_true))**2)

def getClosestSunStations() -> pd.DataFrame:
    '''
    Get the closest sun stations for each location.

    Outputs:
        - cehLocs (pd.DataFrame): DataFrame with closest sun stations
                                 for each location
    '''
    midasLocs = pd.read_csv("data/midas.csv", skiprows=46)
    midasLocs = midasLocs[midasLocs["last_year"] == 2022]
    midasLocs.reset_index(inplace=True)
    cehLocs = pd.read_csv("data/supporting-documents/UKCEHThamesInitiative_SamplingSiteLocations.csv")

    def getLabel(siteName):
        if siteName.split()[0] == "Thames":
            return siteName.split()[-1].lower()
        else:
            return siteName.split()[0].lower()
        
    def closestStation(lat, long):
        loc = (lat, long)
        minDist = 10000000
        minIdx = None
        for idx, row in midasLocs.iterrows():
            dist = hs.haversine(loc, (row["station_latitude"], row["station_longitude"]), unit="km")
            if dist < minDist:
                minDist = dist 
                minIdx = idx 
        return midasLocs.iloc[minIdx]["station_file_name"], minDist

    cehLocs = cehLocs[[c in ["Tm", "Ra","Ch","Oc","TR", "Ke"] for c in cehLocs["Code"]]]
    cehLocs["label"] = [getLabel(s) for s in cehLocs["Site name"]]

    cehLocs["closestStation"] = cehLocs.apply(lambda x : closestStation(x["Latitude"], x["Longitude"])[0], axis=1)
    cehLocs["dist"] = cehLocs.apply(lambda x : closestStation(x["Latitude"], x["Longitude"])[1], axis=1)

    return cehLocs[["label", "closestStation", "dist"]]

class TYPE(Enum):
    CALIBRATION = 1
    VALIDATION = 2
    ALL = 3
