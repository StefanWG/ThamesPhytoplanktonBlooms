import numpy as np
import warnings
from utils import *
from dataProcessing import *
from smoothing import *
from statsmodels.tools.sm_exceptions import InterpolationWarning
from scipy.optimize import curve_fit

warnings.simplefilter('ignore', InterpolationWarning)

PREDICTORS = ["flow", "temp", "sun3d", "silicon", "SRP"]

def decay(x, m, k, d, d2):
    """
    Calculate the decay function value for given parameters.

    Parameters:
    - x: float, independent variable
    - m: float, maximum growth rate
    - k: float, half-saturation constant
    - d: float, death rate
    - d2: float, quadratic death rate

    Returns:
    - float, decay function value
    """
    return m*(x) / (k+x+x**2/d2) - d

def lr(x, a,b,c,d,e):
    """
    Calculate the linear regression function value for given parameters.

    Parameters:
    - x: list of floats, predictor values
    - a, b, c, d, e: floats, regression coefficients

    Returns:
    - float, linear regression function value
    """
    return x[0]*a+x[1]*b+x[2]*c+x[3]*d+x[4]*e

def monodModel(location, tvar, popts=None, predictors=PREDICTORS, split=1, type=TYPE.ALL):
    """
    Monod-based model for predicting phytoplankton growth rates.

    Parameters:
    - location: str, file location of the data
    - tvar: str, target variable for growth prediction
    - popts: dict, optional, pre-calibrated model parameters
    - predictors: list of str, optional, predictor variables
    - split: float, optional, data split ratio for calibration and validation
    - type: TYPE, optional, type of modeling to perform

    Returns:
    - dict or list, predicted growth values depending on the type of modeling
    """
    assert len(predictors) == 5 # This can be modified by must match number of predictors in LR

    data = readData(location)
    smoothed = {v:getInterpolatedData(data["date"], data[v]) for v in predictors}
    split = int(len(smoothed[predictors[0]]) * split)  

    smoothed = {k:smooth(smoothed[k], 2) for k in smoothed.keys()}
    s = getSmoothedData(data["date"], data[tvar], 3)
    growth = np.diff(s) / s[:-1]

    if type == TYPE.CALIBRATION:
        smoothed = {k:smoothed[k][:split] for k in smoothed.keys()}
        growth = growth[:split]

        popts = {"monod":{}, "lr":None}
        preds = []
        for v in smoothed.keys():
            arr = smoothed[v]
            popt, pcov = curve_fit(decay, arr, growth, maxfev=100000,
                        bounds=[(0, 0, 0, 0), 
                                (np.inf,max(arr),2, np.inf)])
            p = decay(arr, *popt)
            popts["monod"][v] = popt
            preds.append(p)

        m=1
        popt, _ = curve_fit(lr, preds, growth, bounds=[(0,0,0,0,0), (m,m,m,m,m)])
        popts["lr"] = popt

        return popts

    elif type == TYPE.VALIDATION:
        assert "monod" in popts 
        assert "lr" in popts
        smoothed = {k:smoothed[k][split:] for k in smoothed.keys()}
        growth = growth[split:]

        monodPreds = []
        for v, popt in popts["monod"].items():
            arr = smoothed[v]
            p = decay(arr, *popt)
            monodPreds.append(p)

        grPreds = lr(monodPreds, *popts["lr"])

        return grPreds

    elif type == TYPE.ALL:
        assert "monod" in popts 
        assert "lr" in popts
        monodPreds = []
        for v, popt in popts["monod"].items():
            arr = smoothed[v]
            p = decay(arr, *popt)
            monodPreds.append(p)

        grPreds = lr(monodPreds, *popts["lr"])

        return grPreds
    