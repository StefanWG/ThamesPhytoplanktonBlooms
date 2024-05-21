from abc import ABC, abstractmethod
from utils import *
from scipy.optimize import curve_fit

# TODO: Add verbose option

class Model(ABC):
    '''
    Abstract class for model. Provides wrapper for 
    model that allows for simple calibration, 
    validation, and plotting. 

    Some things must be defined in order to 
    implement this class:
        - Variables
            - xVariables (list of time series variables included in model)
            - yVariables (list of target variables)
        - Functions
            - model (function that take xData and parameters as input
                     and outputs y data)
            - runModel (function to run model, this can be passed)
    '''
    xData = {}
    yData = {}
    splitWeek = 409 # End of 2018
    sigma = 3 # Choosen based on visual analysis of smoothing
    dates = None
    passYtoModel = False # True if y time series should be passed to the model

    def __init__(self,location, sigma=3, splitWeek = 409):
        self.location = location
        data = readData(location)
        # Initialize self.xData
        for xVar in self.xVariables:
            smoothed = getSmoothedData(data["date"], data[xVar], sigma)
            self.xData[xVar] = smoothed
        # Initialize self.yData
        for yVar in self.yVariables:
            smoothed = getSmoothedData(data["date"], data[yVar], sigma)
            self.yData[yVar] = smoothed

        self.sigma = sigma
        self.splitWeek = splitWeek
        # Intiialize self.dates
        weeksSinceStart = [((i - data["date"][0]).days + 1) // 7 + 1 for i in data["date"]]
        self.dates = [data["date"][0] + pd.Timedelta(weeks=i) for i in range(max(weeksSinceStart))]

        if self.splitWeek > len(data)*0.75:
            self.splitWeek = int(len(data)*0.75)

    # TODO: Pick splitweek smartly - ie. start of year


    @property
    @abstractmethod
    def xVariables(self):
        pass

    @property
    @abstractmethod
    def yVariables(self):
        pass
    
    def calibrate(self):
        '''
        Calibrate model. Uses 0:self.splitWeek for calibration.

        Outputs:
            - calibParameters (dict): Dictionary of yVariables and optimal paramters
        '''
        calibParameters = {}

        xCalib = {x:self.xData[x][:self.splitWeek] for x in self.xVariables}

        for yVar in self.yVariables:
            yCalib = self.yData[yVar][:self.splitWeek]
            if self.passYtoModel:
                xCalib["y"] = yCalib
            popt, pcov = curve_fit(self.model, xCalib,yCalib, maxfev=100000)
            calibParameters[yVar] = popt

        return calibParameters
    
    def validate(self, calibParameters):
        '''
        Validate model. Uses self.splitWeek: for validation.

        Inputs:
            - calibParamters (dict): dict mapping yVariables to optimal paramters,
                                     this is output from self.validate()

        Outputs:
            - validResults (dict): dict mapping yVariables to rmse and Y output
        '''
        validResults = {}
        xValid = {x:self.xData[x][self.splitWeek:] for x in self.xVariables}

        for yVar in self.yVariables:
            yValid = self.yData[yVar][self.splitWeek:]
            if self.passYtoModel:
                xValid["y"] = yValid

            Y = self.model(xValid, *calibParameters[yVar])
            validResults[yVar] = {
                "rmse":rmse(Y, yValid),
                "Y":Y
            }

        return validResults
    
    def plot(self, toCalibrate=True, toValidate=True, output=None):
        '''
        Plot model output for all yVariables. Plot calibration,
        validation, or both.
        '''
        assert toCalibrate or toValidate # One must be true

        calibParameters = self.calibrate()
        xCalib = {x:self.xData[x][:self.splitWeek] for x in self.xVariables}


        if toValidate:
            validResults = self.validate(calibParameters)

        fig = plt.figure(figsize=(8,8))
        years = mdates.YearLocator()

        titles = {
            "diatoms":"Diatoms",  
            "n-chloro":"Nano-chlorophytes", 
            "p-chloro":"Pico-chlorophytes",
            "totCyano":"Cyanobacteria",
            "chl":"Chlorophyll-a"
        }

        fig.tight_layout()
        i = 1
        for yVar in self.yVariables:
            ax = fig.add_subplot(len(self.yVariables), 1, i)
            ax.set_title(titles[yVar])
            i += 1 

            if self.passYtoModel:
                xCalib["y"] = self.yData[yVar][:self.splitWeek]

            if toCalibrate:
                ax.plot(self.dates[:self.splitWeek], self.yData[yVar][:self.splitWeek], "k-", label='Observed')
                ax.plot(self.dates[:self.splitWeek], self.model(xCalib, *calibParameters[yVar]),"r--", label='Model')

            if toValidate:
                ax.plot(self.dates[self.splitWeek:], self.yData[yVar][self.splitWeek:], "k-", label='Observed')
                ax.plot(self.dates[self.splitWeek:], validResults[yVar]["Y"],"r--", label='Model')
            if toCalibrate and toValidate:
                ax.axvline(self.dates[self.splitWeek], color="k", ls=':')
            
            ax.xaxis.set_major_locator(years)
        # Legend
        h, l = ax.get_legend_handles_labels()
        fig.legend(h[:2],l[:2], loc="upper center",ncol=2, bbox_to_anchor=(0.5, .95))
        fig.subplots_adjust(hspace=0.4)

        if output is not None:
            plt.savefig(output)
        else:
            plt.show()

    @abstractmethod
    def model():
        pass

    @abstractmethod
    def runModel():
        pass
