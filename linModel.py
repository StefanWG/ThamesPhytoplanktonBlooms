from model import Model

class LinModel(Model):
    '''
    Implements a model that uses a linear combination 
    of various parameters that is squared (to avoid 
    negative cell counts).

    Model is in the following form:

        $Y = ((\beta•VARS)*z+f)^2$

    where $\beta$ is a vector of coefficients and VARS 
    is the vector of time series data from that time step.
    '''

    xVariables = ["flow","temp","sun3d","SRP","silicon"]
    yVariables = ["diatoms","n-chloro", "p-chloro","totCyano", "chl"]
    
    def model(self, x, a,b,c,d, e,f,z):
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
        F = x["flow"]
        T = x["temp"] # Temp
        S = x["sun3d"] # Sun (3 day avg)
        SRP = x["SRP"] # Soluble reactive phosphorus
        SIL = x["silicon"] # Silicon
        Y = ((a*F + b*T + c*S + d*SRP+e*SIL)*z+ f) **2

        # Y = [max(y,0) for y in Y] include if no exponent
        return Y
    
    def runModel(self):
        xD = {x:self.xData[x] for x in self.xVariables}
        out = self.model(xD, 1,1,1,1,1,1,1)
        

# l = LinModel()
# l.plot()
# l.plot(output="figures/linModel.png")

