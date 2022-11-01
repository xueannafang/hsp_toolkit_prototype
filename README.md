Last update: 01/11/2022

# HSP-toolkits (version 1.0)

* *Solvent Predictor*: Based on the target Hansen solubility parameters (HSPs), propose a list of multi-solvent combination.
* *M Locator*: Predict HSPs of the studied material based on a solubility score.

 <p>
  <img src="https://github.com/xueannafang/hsp-toolkits/blob/main/figs/HSP_toolkit_scheme.png" width=1000>
 </p>
 
## Table of contents
* [General Info](#general-info)
* [*Solvent Predictor*](#solv_pred)
* [*M Locator*](#m_loc)
* [References](#ref)
* [Contribution](#contri)
* [Acknowledgement](#ack)

## General Info

These toolkits are developed with Python 3.7.3 and tested in Windows 10.

The following software are required to be installed:

*Python 3.7, Jupyter Notebook, Microsoft Excel*

The following python packages are required to be installed:

*numpy, pandas, scipy, itertools, abc, os*

## *Solvent Predictor*

Given that HSP of a solvent mixture follows a linear combination of each individual component, researchers can easily calculate the HSP of solvent systems with known components. It is however much more difficult to reverse this process.

The aim of *Solvent Predictor* is to support the solvent suggestion when a desired goal of HSP is known.

The key function is to convert the target HSP into a multi-solvent list based on the requirement of user.

### Set up

- Download the folder of **HSP_SolventPredictor** to local working directory on Windows.
- Open **Solv_pred_class.ipynb** using Jupyter Notebook.

*Please note: **HSP_SolvP.py** contains key calculation process for *Solvent Predictor*. It is imported at the first step. Please be careful to change its name.*
 
### Run *Solvent Predictor*

**Step 1 Import all the related packages:**

```
import numpy as np
import pandas as pd
from scipy.linalg import pinv
import itertools
import abc
import HSP_SolvP as HSP
import os
```

**Step 2 Prepare input spreadsheets using Microsoft Excel:**

- **solvent candidates** (input_solv_sel.xlsx)

 This is the solvent pool for users to choose from.
 
 Solvents appear on this list will be considered as candidates of multi-solvent systems.

 <p>
  <img src=https://github.com/xueannafang/hsp-toolkits/blob/main/figs/input_solv_sel_example.png>
  <em>Example of solvent candidates list<em>
 </p> 


   Users can edit this form depending on their requirement.
   
   *To remove an entry*:
   - Select the whole row and delete.
   - Update the number in "Group" column.
   **Note: The "Group" number must be continious integers starting from 1 to the total number of solvents.**
   **Empty, discontinous, or overflowed numbers in "Group" column will lead to error.**


**Step 3**
Specify calculation parameters:

**n**
- maximum number of solvents in each multi-solvent combination
- default = 2
- must be an integer

**target HSP**
- contains three parameters in the order of delta_D, delta_P, delta_H
- must be float



- repeated calculation times for purturbation applied on the target HSP matrix


Upload solvent candidate list (input_solv_sel.xlsx) and database (db.xlsx)

```
class SolvPred():
    
    def __init__(self, input_solv, db):
        self.pred = HSP.SolvPredictor(input_solv, db)

    def mix_pred(self, n = 2, rep_time = 50, std = 0.1, tol_pred = 1, red_tol = 0.01):
        self.pred.run_all(n, 18.0, 1.4, 2.0, rep_time = rep_time, std = std, tol = tol_pred, red_tol = red_tol)

sp = SolvPred(r'input_solv_sel.xlsx', r'db.xlsx')

sp.mix_pred()

```


## Contribution
This project is developed by
- Xue Fang (School of chemistry, University of Bristol, UK)
- Bo Gao (School of physics, University of Bristol, UK)

with guidance from
- Prof. Charl FJ Faul (School of chemistry, University of Bristol, UK)
- Dr Natalie Fey (School of chemistry, University of Bristol, UK)
- Dr Ella Gale (School of chemistry, University of Bristol, UK)

Experimental test data is contributed by the CMP team of the Faul research group.

## Acknowledgements

University of Bristol

Chinese Scholarship Council

Royal Society of Chemistry
