# River Thames Algae Blooms - Dissertation
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

### Model 

The model is given by the following equation:

$$Y_t=(\beta_6 * (\boldsymbol{\beta} \cdot \bf{V}) + \beta_7)^2 $$

where $\boldsymbol{\beta} = \\{ \beta_1, \beta_2,\beta_3, \beta_4,\beta_5 \\}$
is a vector of coefficients and $\bf{V}$ is a vector of temperature,
flow, sun, phosphorus, and silicon time series data at time $t$.

<img src="figures/model.png" alt="model" width="600"/>


