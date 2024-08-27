from abc import ABC, abstractmethod
import datetime
from dataProcessing import readData
from utils import *
import matplotlib.pyplot as plt 
import matplotlib.dates as mdates
import kneed 
from smoothing import *
from scipy.signal import find_peaks, peak_widths
import seaborn as sns


def getCalibrationBlooms():
    """
    Get the calibration blooms by replicating blooms from Bowes et al. (2024).

    Returns:
        numpy.ndarray: The calibration blooms.
    """
    data = readData("runnymede")
    data = data[data['date'] <=datetime.datetime(2018, 12, 31)]
    flowSmoothed = getSmoothedData(data["date"], data["flow"], SIGMA)
    tempSmoothed = getSmoothedData(data["date"], data["temp"], SIGMA)
    sun3dSmoothed = getSmoothedData(data["date"], data["sun3d"], SIGMA)
    diatomsSmoothed = getSmoothedData(data["date"], data["diatoms"], 1)
    tempThresh = np.logical_and(tempSmoothed >= 10.1, tempSmoothed <= 23)
    flowThresh = np.logical_and(flowSmoothed >= 0.7, flowSmoothed <= 100)
    sunThresh = sun3dSmoothed >=15000/24
    thres = np.logical_and(np.logical_and(tempThresh, flowThresh), sunThresh)
    thres = diatomsSmoothed >= np.nanpercentile(diatomsSmoothed, 75)
    return thres

class Blooms(ABC):
    """
    Abstract base class for bloom detection algorithms.
    """

    end_date = datetime.datetime(2018, 12, 31)
    blooms = None
    starts = None
    ends = None 
    peaks = None

    def __init__(self, location, tvar):
        """
        Initialize the Blooms object.

        Args:
            location (str): The location of the data.
            tvar (str): The variable to use for bloom detection.
        """
        self.location = location
        self.tvar = tvar

        self.data = readData(self.location)
        self.data = self.data[self.data['date'] <= self.end_date]

        self.smoothed = getSmoothedData(self.data["date"], self.data[tvar], SIGMA)

    @abstractmethod
    def getBlooms(self):
        """
        Abstract method to get the bloom detection results.
        """
        pass
    
    @abstractmethod
    def calibrate(self):
        """
        Abstract method to calibrate the bloom detection algorithm.
        """
        pass

    def getThresholds(self, confInterval=100, vars=["temp", "flow", "sun3d"]):
        """
        Get the thresholds for the specified variables.

        Args:
            confInterval (int, optional): The confidence interval for the thresholds. Defaults to 100.
            vars (list, optional): The variables to get thresholds for. Defaults to ["temp", "flow", "sun3d"].

        Returns:
            dict: A dictionary containing the thresholds for each variable.
        """
        assert self.blooms is not None
        self.thresholds = {}
 
        for var in vars:
            weeksSinceStart = getWeeksSinceStart(self.data["date"])
            # Create array to extract thresholds from
            varMap = {a:b for a,b, in zip(weeksSinceStart, self.data[var])}    
            v = np.array([varMap[i+1] if i+1 in varMap else np.nan for i in range(max(weeksSinceStart))])
            idxs = np.where(self.blooms)[0]
            self.thresholds[var] = {
                "min":np.nanpercentile(v[idxs], (100-confInterval) / 2), 
                "max":np.nanpercentile(v[idxs], (100-confInterval) / 2 + confInterval)
            }

        return self.thresholds

    def validate(self, location, tvar):
        """
        Validate the bloom detection algorithm using a different location and variable.

        Args:
            location (str): The location of the validation data.
            tvar (str): The variable to use for validation.

        Returns:
            Blooms: A new Blooms object with the validation data.
        """
        valid = self.__class__(location, tvar)
        return valid

class KneeIdentification(Blooms):     
    """
    Class for bloom detection using the knee identification method.
    """

    knee = None
    optSigma = None

    def __init__(self, location, tvar):
        """
        Initialize the KneeIdentification object.

        Args:
            location (str): The location of the data.
            tvar (str): The variable to use for bloom detection.
        """
        super().__init__(location, tvar)

    def getBlooms(self, sigma=None, S=1):
        """
        Get the bloom detection results using the knee identification method.

        Args:
            sigma (int, optional): The sigma value for smoothing. Defaults to None.
            S (int, optional): The sensitivity parameter. Defaults to 1.

        Returns:
            numpy.ndarray: The bloom detection results.
        """
        if sigma is None:
            if self.optSigma is not None:
                sigma = self.optSigma 
            else:
                sigma = 10

        knee = self.getKneedleThreshold(sigma=sigma, S=S)
        # ...
        # Rest of the code
        # ...
        self.blooms = self.smoothed >= knee
        self.starts = []
        self.ends = []
        prev = False
        for i, v in enumerate(self.blooms):
            if not prev and v:
                self.starts.append(i)
            elif prev and not v:
                self.ends.append(i)
            prev = v
        if prev: 
            self.ends.append(len(self.smoothed)-1)
        self.peaks = [(ma+mi)//2 for ma, mi in zip(self.starts, self.ends)]
   
        return self.blooms

    def getKneedleThreshold(self,sigma=5,S=1):
        sorted = np.sort(self.smoothed)
        sorted = sorted[~np.isnan(sorted)]
        s_sorted = smooth(sorted, sigma)
        perc = np.linspace(0, 1, len(s_sorted))
        kneedle = kneed.KneeLocator(s_sorted, perc, curve="concave", S=S)
        self.knee = kneedle.knee
        return kneedle.knee 
    
    def calibrate(self):
            """
            Calibrates the model by finding the optimal sigma value based on the F-statistic.

            Returns:
                float: The optimal sigma value for the calibration.
            """
            calib = getCalibrationBlooms(self.location, self.tvar)
            def kneeCalib(data, sigma):
                thresh = self.getKneedleThreshold(sigma=sigma)
                modelBlooms = self.smoothed >= thresh
                f = fStat(calib, modelBlooms)
                return f

            # Maximize F stat
            sigmas = np.arange(1,50)
            res = [kneeCalib(self.data, i) for i in sigmas]
            idx = np.where(res==max(res))[0][0]
            # print(f"F-Stat inflection: {round(max(res), 3)}, Sigma: {sigmas[idx]}")
            self.optSigma = sigmas[idx]
            return self.optSigma

class PeakIdentification(Blooms):
    """
    Class for identifying peaks and calibrating bloom detection.

    Args:
        location (str): The location of the data.
        tvar (str): The time variable.

    Attributes:
        optProm (int or None): The optimal prominence value for peak detection.

    """

    optProm = None

    def __init__(self, location, tvar):
        super().__init__(location, tvar)

    def getBlooms(self, prom=None):
        """
        Get the bloom detection results.

        Args:
            prom (float or None, optional): The prominence value for peak detection. If None, the optimal prominence value will be used. Defaults to None.

        Returns:
            list: A list of boolean values indicating the presence of blooms.

        """
        if prom is None:
            if self.optProm is not None:
                prom = self.optProm
            else:
                prom = 10
        self.peaks, _ = find_peaks(self.smoothed, prominence=np.max(self.smoothed) / prom)
        widths = peak_widths(self.smoothed, self.peaks)
        self.starts = np.round(widths[2]).astype(int)
        self.ends = np.round(widths[3]).astype(int)

        self.blooms = [False for _ in range(len(self.smoothed))]
        for i in range(len(self.peaks)):
            for x in range(self.starts[i], self.ends[i] + 1):
                self.blooms[x] = True
        return self.blooms

    def calibrate(self):
        """
        Calibrate the peak detection by maximizing the F-statistic.

        Returns:
            int: The optimal prominence value for peak detection.

        """
        calib = getCalibrationBlooms(self.location, self.tvar)

        def peakCalib(prom):
            modelBlooms = self.getBlooms(prom=prom)
            f = fStat(calib, modelBlooms)
            return f

        # Maximize F stat
        proms = np.arange(1, 50)
        res = [peakCalib(i) for i in proms]
        idx = np.where(res == max(res))[0][0]
        # print(f"F-Stat peaks: {max(res)}, Prominence: {proms[idx]}")
        self.optProm = proms[idx]
        return self.optProm

