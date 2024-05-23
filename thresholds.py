
import numpy as np
from smoothing import *
import scipy
import kneed

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

def getKneedleThreshold(data, S = 1):
    s = np.sort(data)
    s = s[~np.isnan(s)]
    smoothed = smooth(s, 5)
    perc = np.linspace(0, 1, len(smoothed))
    kneedle = kneed.KneeLocator(smoothed, perc, curve="concave", S=S)
    return kneedle.knee

def getInflectionThreshold(data):
    s = np.sort(data)
    s = smooth(s, 50) # NOTE: This parameter is super important
    deriv = np.diff(s)
    idx = np.where(deriv == max(deriv))[0][0]
    return s[idx]

def getBlooms(dates, vec, method="find_peaks", sigma=None, type="gaussian", kneedleS = 1):
    '''
    Get dates during with blooms and cessations occured.
    This is given by:
        - peaks using the scipy find_peaks algorithm.
        - the knee of the PDF using Kneedle algorithm

    Inputs:
        - dates (array): list of dates corresponding to time series data
        - vec (array): time series data
        - method (str): method of bloom finding; either "find_peaks" or
                        "kneedle" (default: "find_peaks)
        - sigma (int): sigma for gaussian curve (default: None)
        - type (str): type of filter to use for convolution; either
                      "gaussian" or "uniform" (default: "gaussian)
        - kneedleS (int): S parameter for kneed.KneeLocator (default: 1)

    Outputs:
        - smoothed (array): array of smoothed and gap filled time series data
        - peaks (array): array of peaks (this corresponds to where bloom becomes cessation)
        - starts (array): array of bloom starting points
        - ends (array): array of cessation ending points


    '''
    assert method in ["find_peaks", "kneedle", "inflection"]
    smoothed = getSmoothedData(dates, vec, sigma, type)

    if method == "find_peaks":
        peaks, _ = scipy.signal.find_peaks(smoothed, prominence=np.max(smoothed) /20)
        widths = scipy.signal.peak_widths(smoothed, peaks)
        starts = np.round(widths[2]).astype(int)
        ends = np.round(widths[3]).astype(int)
    elif method == "kneedle":
        knee = getKneedleThreshold(smoothed, kneedleS)
        blooms = smoothed >= knee
        starts = []
        ends = []
        prev = False
        for i, v in enumerate(blooms):
            if not prev and v:
                starts.append(i)
            elif prev and not v:
                ends.append(i)
            prev = v
        peaks = [(ma+mi)/2 for ma, mi in zip(starts, ends)]
        if len(ends) < len(starts):
            ends.append(len(smoothed)-1)
    elif method == "inflection":
        thresh = getInflectionThreshold(smoothed)
        blooms = smoothed >= thresh
        starts = []
        ends = []
        prev = False
        for i, v in enumerate(blooms):
            if not prev and v:
                starts.append(i)
            elif prev and not v:
                ends.append(i)
            prev = v
        peaks = [(ma+mi)/2 for ma, mi in zip(starts, ends)]
        if len(ends) < len(starts):
            ends.append(len(smoothed)-1)


    return smoothed, peaks, starts, ends

