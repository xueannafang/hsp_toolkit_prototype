Last update on 07/11/2022

# M Locator

As indicated by its name, *M Locator* locates "*M*" in the Hansen space. It is a top-down approach aiming at obtaining the HSP of as-studied material (M) using a solubility score (*w<sub>i</sub>*) measured from a series of solvents.

*M Locator* can be useful for complicated systems whose exact chemical component is hard to confirm. This could be, for example, polymeric systems with high dispersity, mixtures of isomers or analogous, etc., where bottom-up method like group contribution may not be easy to be carried out. 

<p>
  <img src=https://github.com/xueannafang/hsp-toolkits/blob/main/figs/m_loc_sch.png width=500>
  </p>

Part of the *M Locator* is connected to [*Solvent Predictor*](https://github.com/xueannafang/hsp-toolkits/edit/main/HSP_SolventPredictor/solv_pred_readme.md) to provide multi-solvent suggestion based on the as-predicted HSP of the material.

The solubility indicator applied in this version is based on UV-Vis spectroscoy features.

The studied materials in this case are based on [conjugated microporous polymers (CMP)](https://www.sciencedirect.com/topics/engineering/conjugated-microporous-polymer).


## Before start

[Python](https://www.python.org/),
[Jupyter notebook](https://jupyter.org/)
and [Microsoft Excel](https://www.microsoft.com/en-us/microsoft-365/excel) are required to be installed before using this tool.

[numpy](https://numpy.org/),
[pandas](https://pandas.pydata.org/), 
[scipy](https://docs.scipy.org/doc/), 
[mpl_toolkits](https://matplotlib.org/1.3.0/mpl_toolkits/index.html), and 
[matplotlib.pyplot](https://matplotlib.org/3.5.3/api/_as_gen/matplotlib.pyplot.html) are required to be installed in python.

**OS:** Windows 10

## Set up

- Download the folder of [**HSP_MLocator**](https://github.com/xueannafang/hsp-toolkits/tree/main/HSP_MLocator) to local working directory on Windows.
- Open **M_loc_class.ipynb** using Jupyter Notebook.

## Run *M Locator*

### Step 1 Load all the related packages:

Execute the first cell by **shift + enter**
```
import numpy as np
import pandas as pd
from scipy.linalg import pinv
import itertools
import HSP_SolvP as HSP
import HSP_M_loc as LOCM
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import os
```
### Step 2 Prepare input spreadsheets using Microsoft Excel:

Note: Explanation of how to use *input_solv_sel.xlsx* and *db.xlsx* is described in **Step 2** of [*Solvent Predictor*](https://github.com/xueannafang/hsp-toolkits/edit/main/HSP_SolventPredictor/solv_pred_readme.md).

- Solubility scores (*input_mloc_data.xlsx*)

This spreadsheet contains experimental-measured results that serve as the solubility score when optimising the location of M.

Users need to update this spreadsheet based on their own experimental set up.

Hereby, the column of "Group", "No_db", "CAS", "Solvent" follow the same method as described in [*Solvent Predictor*](https://github.com/xueannafang/hsp-toolkits/edit/main/HSP_SolventPredictor/solv_pred_readme.md).

In this version, all the "ratio" are set as 100%, suggesting the neat solvent has been applied.

(The future versions will start to include multi-solvent input, where both the "Group" and "ratio" will work together, but for now, don't worry about that.)






