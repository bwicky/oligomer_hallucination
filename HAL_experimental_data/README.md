# Experimental data

Experimental data collected on designs reported in the [oligomer hallucination paper](200~https://www.biorxiv.org/content/10.1101/2022.06.09.493773v1.full).

`HAL_exp_data.csv` contains ...

`HAL_exp_data.h5` contains ... as well as all SEC and CD chromatograms. In order to load that dataframe, do:

```
import pandas as pd
df = pd.read_hdf('HAL_exp_data.h5', 'df')
```

Definitions
