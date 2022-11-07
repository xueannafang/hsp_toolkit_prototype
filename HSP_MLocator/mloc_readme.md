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

(In the future versions, there will be functions of multi-solvent input, where both the "Group" and "ratio" will work together, but for now, don't worry about that.)

The "Indicator" column is the solubility indicator, which in this case is maximum characteristic UV-Vis absorbance data.

*Please make sure all the data filled in this column must be non-negative float.*

### Step 3 Specify arguments in run_all()

*run_all()* contains 9 arguments, among which **alpha, n_max, tol_mop** belongs to the optimisation of HSP of M, and **n, rep_time, std, tol_pred, red_tol** belongs to *Solvent Predictor* (which will not be introduced in this document).

```
run_all(self, alpha = 0.001, n_max = 1000000, tol_mop = 0.0001, n = 3, rep_time = 50, std = 0.1, tol_pred = 0.1, red_tol = 0.01)
```

Given that less Hansen distance suggests better solvent, which should ideally have higher solubility indicator, the key optimisation process is based on minimise the overall weighted Hansen distance between M and all the tested solvents, where the weight of Hansen distance is the solubility score we mentioned before.

We use [gradient descent](https://ml-cheatsheet.readthedocs.io/en/latest/gradient_descent.html) to find the optimised point.

**alpha**

- learning rate
- default = 0.001
- too big may cause the fail of converge
- too small may increase the calculation time and interation times

**n_max**
- maxmium iteration times
- default = 1000000
- too small may cause fail of converge

**tol_mop**
- an extreme small number served as the threshold below which the system is regarded as converged
- default = 0.0001

## Step 4 Upload input files and run

Upload *input_mloc_data.xlsx*, *db.xlsx* and *input_solv_sel.xlsx* into the last second line in the cell below.

Execute this cell.

```
class MLocSolvPred():

    def __init__(self, exp_result, input_solv, db):
        self.mop = LOCM.MLoc(exp_result, db)
        print(self.mop.folder_path)
        self.pred = HSP.SolvPredictor(input_solv, db, folder_name = self.mop.folder_path)
    
    def run_all(self, alpha = 0.001, n_max = 1000000, tol_mop = 0.0001, n = 3, rep_time = 50, std = 0.1, tol_pred = 0.1, red_tol = 0.01):
        cord_result = self.mop.run_all(alpha = alpha, n_max = n_max, tol = tol_mop)
        self.pred.run_all(n, cord_result[0][0], cord_result[0][1], cord_result[0][2], rep_time = rep_time, std = std, tol = tol_pred, red_tol = red_tol)

mp = MLocSolvPred(r'input_mloc_data.xlsx', r'input_solv_sel.xlsx', r'db.xlsx')
mp.run_all()
```

## Output examples

Once the calculation has finished, you will see a rough calculation report in the notebook:
<p>
  <img src=https://github.com/xueannafang/hsp-toolkits/blob/main/figs/op_mloc_notebook.png>
  </p>

where the folder path, initial guess, iteration steps as well as a 3D plot containing M and all the applied solvents will be presented.

At the mean time, a folder named *mloc_data* is created.

(Note that if the *input_mloc_data.xlsx* has been renamed, this folder name would be anything between *input_* and *.xlsx*) 

The output folder contains

- a png file of the 3D scheme shown in the notebook
- a spreadsheet with all the HSP involved in the calculation and comparison of solvent quality using Hansen distance
- a log file of calculation details and optimisation results
- the full output based on *Solvent Predictor*, where the target is the HSP of M

<p>
  <img src=https://github.com/xueannafang/hsp-toolkits/blob/main/figs/op_mloc_files.png>
  </p>

**\*HSP_compare.xlsx**

This is the most important output that the user need to read.

<p>
  <img src=https://github.com/xueannafang/hsp-toolkits/blob/main/figs/op_mloc_compare.png>
  </p>
 
 In the "MLoc" row, the HSP of material is list there, with a suggestion of best solvent whose Hansen distance from M is the smallest.
 
 The rest of solvents are also summarised in this spreadsheet, with their HSP list in the **D, P, H, T** columns.
 
 The difference between each HSP sub-parameter is calculated and list in the **dD, dP, dH, dT** columns.
 
 (This design is to make it flexible for user to choose the most interested sub-parameter to understand the contribution of different molecular interaction.)
 
 The final column is the calculated Hansen distance (**R**). Solvents are ranked in the ascending sequence of **R**.
 
 Hereby we keep the Hildebrand (or total) solubility parameter (**T**), which is another useful solubility parameter but does not separate the three interactions.
