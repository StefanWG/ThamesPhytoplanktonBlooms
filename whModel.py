from model import Model
import numpy as np

class WHModel(Model):
    '''
    Implements algorithm from Whitehead and Hornberger (1984).
    Original algorithm uses chlorophyll-a instead of cell counts.
    '''

    #TODO: Why are there overflow issues?

    xVariables = ["flow", "sun3d"]
    yVariables = ["diatoms","n-chloro", "p-chloro","totCyano"]
    passYtoModel = True

    def getdydt(self, Q, cc, sr, k2, k3, k4, k5, k6, k7):
        '''
        Original algorithm is commented. Removed first two terms 
        because we are predicting cell counts at a specific location
        (i.e. an instantaneous solution) and Qd=Qu and chl_d=chl_u.

        Inputs:
            - Q (float): discharge (flow rate)
            - cc (float): cell count at previous time step
            - sr (float): solar radition
            - k2 (float): tunable parameter
            - k3 (float): tunable parameter
            - k4 (float): tunable parameter
            - k5 (float): tunable parameter
            - k6 (float): tunable parameter
            - k7 (float): tunable parameter

        Outputs:
            - dydt (float): change in cell count
        '''
        #TODO: More detail in doc string
        #TODO: Can we add water temperature to the model?

        # return k1*Qu*chl_u \
        # - k1*Qd*chl_d \
        # - k2*chl_d \
        # + k3*sr / Qd \
        # * (k4 / (k4 + Qd**k5)) \
        # * (sr / k6)**k7 \
        # * np.exp(1 - (sr / k6)**k7)
        return - k2*cc \
        + k3*sr / Q \
        * (k4 / (k4 + Q**k5)) \
        * (sr / k6)**k7 \
        * np.exp(1 - (sr / k6)**k7)

    def model(self, x, k2, k3, k4, k5 , k6, k7):
        '''
        Implementation of Whitehead and Hornberger (1984)
        model. 

        Inputs:
            - x (dict): dict mapping xVariables to time series data
            - k2 (float): tunable parameter
            - k3 (float): tunable parameter
            - k4 (float): tunable parameter
            - k5 (float): tunable parameter
            - k6 (float): tunable parameter
            - k7 (float): tunable parameter

        Outputs:
            - Y (np.array): cell count time series
        '''
        F = x["flow"]
        S = x["sun3d"]
        yData = x["y"]

        Y = [yData[0]]
        for i in range(len(yData)-1):
            dydt = self.getdydt(F[i], Y[-1], S[i], k2, k3, k4, k5, k6, k7)
            Yn = Y[-1] + dydt
            Y.append(Yn)

        return np.array(Y)

    def runModel(self):
        pass

wh = WHModel()
wh.plot(output="figures/whModel.jpg")