# River Thames Phytoplankton Blooms - Dissertation
## Stefan Walzer-Goldfeld

### Multi-Parameter Plot
<img src="figures/multiParamCells.jpg" alt="drawing" width="600"/>

### Time Series Plots
**Diatoms**

<img src="figures/diatoms.jpg" alt="diatoms" width="600"/>

**Nano-Chlorophytes**

<img src="figures/nchloro.jpg" alt="nchloro" width="600"/>

**Pico-Chlorophytes**

<img src="figures/pchloro.jpg" alt="pchloro" width="600"/>

**Total Cyanobacteria**

<img src="figures/totCyano.jpg" alt="totCyana" width="600"/>

**Total Phytoplankton**

<img src="figures/totPhyto.jpg" alt="totPhyto" width="600"/>

### Simple Model 

This model is given by the following equation:

$$Y_t=(\beta_6 * (\boldsymbol{\beta} \cdot \bf{V}) + \beta_7)^2 $$

where $\boldsymbol{\beta} = \\{ \beta_1, \beta_2,\beta_3, \beta_4,\beta_5 \\}$
is a vector of coefficients and $\bf{V}$ is a vector of temperature,
flow, sun, phosphorus, and silicon time series data at time $t$.

<img src="figures/linModel.png" alt="linModel" width="600"/>

### Whitehead and Hornberger (1984) Model 

This model is an implementation of the model developed by Whitehead and Hornberger (1984).
The original model uses chlorophyll-a as a proxy for phytoplankton cell counts, however 
this does not. The original model also uses upstream and downstream data, while this 
implementation uses data from a single location for both upstream and downstream points 
(i.e. we are computing an instantaneous solution).

<img src="figures/whModel.jpg" alt="whModel" width="600"/>

Note: in some cases, the optimization algorithm does not converge; this is a
work in progress.


