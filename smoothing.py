
import numpy as np

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

    If gap is large than the size of the filter (6*sigma),
    then gap filling will give values of 0.

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
    return np.array(smoothed)

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
    toSmooth = np.array([vMap[i+1] if i+1 in vMap else 0 for i in range(max(weeksSinceStart))], dtype=float)
    toSmooth[np.isnan(toSmooth)] = 0
    # Smooth
    smoothed = np.array(smooth(toSmooth, sigma, type))
    return smoothed

