import numpy as np
import warnings
import pandas as pd
from utils import *
from dataProcessing import *
from smoothing import *

from statsmodels.tools.sm_exceptions import InterpolationWarning
from dataProcessing import *
from statsmodels.regression.linear_model import OLS
from sklearn.metrics import root_mean_squared_error

import kneed


warnings.simplefilter('ignore', InterpolationWarning)

def lagModel(location, tvar, popts=None, split=1, type=TYPE.ALL, start=1, end=2):
    """
    Generate lagged model predictions based on the provided data.

    Parameters:
    - location (str): The location of the data.
    - tvar (str): The variable to be used for modeling.
    - popts (list, optional): The list of model parameters. Defaults to None.
    - split (float, optional): The ratio of data to be used for training. Defaults to 1.
    - type (TYPE, optional): The type of prediction to generate. Defaults to TYPE.ALL.
    - start (int, optional): The starting lag value. Defaults to 1.
    - end (int, optional): The ending lag value. Defaults to 2.

    Returns:
    - preds (array-like): The predicted values based on the lagged model.

    """
    # Note: predictions begin after week 'end'
    if popts is not None:
        end = start - 1 + len(popts)

    def getDF(start, end):
        """
        Generate a lagged dataframe based on the provided start and end lag values.

        Parameters:
        - start (int): The starting lag value.
        - end (int): The ending lag value.

        Returns:
        - df (DataFrame): The lagged dataframe.

        """
        cellConcs = getSmoothedData(data["date"], data[tvar], 3)
        gr = np.diff(cellConcs) / cellConcs[:-1]
        df = pd.DataFrame({"gr": gr})
        for i in range(start, end + 1):
            df[f"lag{i}"] = df["gr"].shift(i)
        df.dropna(inplace=True)
        df = df.reset_index().drop(columns=["index"])
        return df

    data = readData(location)
    df = getDF(start, end)

    split = int(len(df) * split)

    if type == TYPE.ALL:
        preds = (df.drop(columns=["gr"]) * np.array(popts)).sum(axis=1)
        return preds
    elif type == TYPE.VALIDATION:
        df = df[split:]
        preds = (df.drop(columns=["gr"]) * np.array(popts)).sum(axis=1)
        return preds
    elif type == TYPE.CALIBRATION:
        errs = []
        for lagEnd in range(start, 9):
            df = getDF(start, lagEnd)
            train = df[:split]

            model = OLS(train["gr"], train.drop(columns=["gr"]))
            res = model.fit()
            preds = res.predict(train.drop(columns=["gr"]))
            err = root_mean_squared_error(preds, train["gr"])
            errs.append(err)

        k = kneed.KneeLocator(np.arange(start, 9), errs, curve="convex", direction="decreasing")
        knee = k.knee

        df = getDF(start, knee + 1)
        df = df.reset_index().drop(columns=["index"])

        train = df[:split]

        model = OLS(train["gr"], train.drop(columns=["gr"]))
        res = model.fit()

        return res.params
