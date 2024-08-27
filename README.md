# River Thames Phytoplankton Blooms - Dissertation
## Stefan Walzer-Goldfeld

This repository contains the follwoing tools:
- Methods for Identifying Blooms From Time Series Data 
    - Knee of CDF Method
    - Peak Identification Method
- Predictive Lag Model
- Mechanistic Monod-based model

Below are examples of how to use each tool. 

### Bloom Identification and Thresholds
The knee method is used in the follow manner: 

```
from blooms import KneeIdentification

location = "runnymede"
target = "diatoms"
kneeModel = KneeIdentifcation(location, target)
kneeModel.calibrate()
blooms = kneeModel.getBlooms() # np.array; 1 if week is a bloom, 0 otherwise
thresholds = kneeModel.getThresholds() # dict
```

The peak identificaiton method uses the `PeakIdentification`
class and is used in the same manner.

### Lag Model
The lag model is used in the follow manner:

```
from lagModel import lagModel 
from utils import TYPE

start = 1
location = "runnymede"
target = "diatoms"
split = 0.7
popts = lagModel(location, target, start=start, split=split,
                 type=TYPE.CALIBRATION)

preds = lagModel(location, target, start=start, split=split,
                 type=TYPE.VALIDATION) # Predictions
```

### Monod-based Model
```
from mondModel import mondModel 
from utils import TYPE

location = "runnymede"
target = "diatoms"
split = 0.7
popts = lagModel(location, target, split=split,
                 type=TYPE.CALIBRATION)

preds = lagModel(location, target, split=split,
                 type=TYPE.VALIDATION) # Predictions
```

### Further Information 
All tools assume that data can be read using the read
using the `readData` function. This requires that 
data has been processed into a specific format. 
This is done from UKCEH files using 
`python3 driver.py -c`.

All locations are found in `utils.LOCATIONS` 
and phytoplankton groups in `utils.TARGETVARS`.
For any help with the provided code or with plotting,
feel free to raise an issue.



