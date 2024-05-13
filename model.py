import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
from copy import deepcopy
from utils import *
from scipy.optimize import curve_fit

plt.rcParams['text.usetex'] = True

import warnings
warnings.filterwarnings("ignore")

sigma = 3 # Choosen with visual analysis of curve

def model(x, a,b,c,d, e,f,z):
    '''
    Simple model in the form of $Y = ((\beta•VARS)*z+f)^2$
    where $\beta$ is a vector of coefficients and VARS 
    is the vector of time series data from that time step.

    Notes: 
        - Model has no dependence on current or initial cell counts.
        - Square in model ensure that cell count cannot go below 0.

    Inputs:
        - x (np.array): 2d array of time series data (shape: (5,n))
            - x[0]: flow data
            - x[1]: temperature data
            - x[2]: sun data (3-day average)
            - x[3]: soluble reactive phosphorus data
            - x[4]: silicon data
        - a (float): flow coefficient
        - b (float): temp coefficient
        - c (float): sun coefficient
        - d (float): srp coefficient
        - e (float): silicon coefficient
        - f (float): additive value
        - z (float): magnitude coefficient

    Outputs:
        - Y (array): array of modelled cell counts
    '''
    F = x[0] # Flow
    T = x[1] # Temp
    S = x[2] # Sun (3 day avg)
    SRP = x[3] # Soluble reactive phosphorus
    SIL = x[4] # Silitcon
    Y = ((a*F + b*T + c*S + d*SRP+e*SIL)*z+ f)**2
    return Y
    

# Import data
data = readData()
flowSmoothed = getSmoothedData(data["date"], data["flow"], sigma)
tempSmoothed = getSmoothedData(data["date"], data["temp"], sigma)
sun3dSmoothed = getSmoothedData(data["date"], data["sun3d"], sigma)
srpSmoothed = getSmoothedData(data["date"], data["SRP"], sigma)
siliconSmoothed = getSmoothedData(data["date"], data["silicon"], sigma)

weeksSinceStart = [((i - data["date"][0]).days + 1) // 7 + 1 for i in data["date"]]
dates = [data["date"][0] + pd.Timedelta(weeks=i) for i in range(max(weeksSinceStart))]
splitWeek = 409 # End of 2018


fig = plt.figure(figsize=(6,8))
fig.tight_layout()
i = 1
# Titles for plotting
titles = {
    "diatoms":"Diatoms",  
    "n-chloro":"Nano-chlorophytes", 
    "p-chloro":"Pico-chlorophytes",
    "totCyano":"Cyanobacteria",
}
# Initial values - picked ones such that curve_fit finds a solution
p0 = {
    "diatoms":[-1,1,1,1,1,1,1],  
    "n-chloro":[-11,1,1,1,1,1,1], 
    "p-chloro":[-11,1,1,1,1,1,1],
    "totCyano":[-15,1,1,1,1,1,1]
}

for y in ["diatoms", "p-chloro", "n-chloro", "totCyano"]:
    ax = fig.add_subplot(4,1,i)
    i += 1
    # Set up x and y data
    ySmoothed = getSmoothedData(data["date"], data[y], sigma)
    xdata = np.array([[a,b,c,d,e] for a,b,c,d,e in zip(flowSmoothed, tempSmoothed, sun3dSmoothed,   
                                                       srpSmoothed,siliconSmoothed)]).transpose()
    xCalib, xValid = xdata[:,:splitWeek], xdata[:,splitWeek:]
    yCalib, yValid = ySmoothed[:splitWeek], ySmoothed[splitWeek:]

    # Calibrate
    popt, pcov = curve_fit(model, xCalib,yCalib, p0=p0[y], maxfev=100000)

    ax.plot(dates[:splitWeek], yCalib, "k-", label='Observed')
    ax.plot(dates[:splitWeek], model(xCalib, *popt),"r--", label='Model')

    # Validate and plot
    ax.set_title(titles[y])
    ax.plot(dates[splitWeek:], yValid, "k-", label='Observed')
    ax.plot(dates[splitWeek:], model(xValid, *popt), "r--", label='Mod - V')
    ax.axvline(dates[splitWeek], color="k", ls=':')

# Create legend
h, l = ax.get_legend_handles_labels()
fig.legend(h[:2],l[:2], loc="upper center",ncol=2, bbox_to_anchor=(0.5, .95))
fig.subplots_adjust(hspace=0.3)

plt.savefig("figures/model.png")
plt.show()
