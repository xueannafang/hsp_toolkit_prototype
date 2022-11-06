Last update on 06/11/2022

# Solvent Predictor


Given that HSP of a solvent mixture follows a linear combination of each individual component, researchers can easily calculate the HSP of solvent systems with known components. It is however much more difficult to reverse this process.

<p>
 <img src="https://github.com/xueannafang/hsp-toolkits/blob/main/figs/bg_solvpred.png" width=1000>
 </p>

The aim of *Solvent Predictor* is to support the solvent suggestion when a desired goal of HSP is known.

The key function is to convert the target HSP into a multi-solvent list based on the requirement of user.

## Before start

This is a beginner-friendly version that does not require too much coding background.

However, to avoid troubles as possible, please be extra careful with
- Removing or changing file names of .py documents;
- Using space or other special characters in input spreadsheets, especially in "CAS No.";
- Editing title line of input spreadsheets.

[Python](https://www.python.org/),
[Jupyter notebook](https://jupyter.org/)
and [Microsoft Excel](https://www.microsoft.com/en-us/microsoft-365/excel) are required to be installed before using this tool.

[numpy](https://numpy.org/),
[pandas](https://pandas.pydata.org/), and
[scipy](https://docs.scipy.org/doc/) are required to be installed in python.


## Set up

- Download the folder of [**HSP_SolventPredictor**](https://github.com/xueannafang/hsp-toolkits/tree/main/HSP_SolventPredictor) to local working directory on Windows.
- Open **Solv_pred_class.ipynb** using Jupyter Notebook.

## Run *Solvent Predictor*

### Step 1 Load all the related packages:

Execute the first cell by **shift + enter**

```
import numpy as np
import pandas as pd
from scipy.linalg import pinv
import itertools
import abc
import HSP_SolvP as HSP
import os
```

### Step 2 Prepare input spreadsheets using Microsoft Excel:

- **solvent candidates** (input_solv_sel.xlsx)

 This is the solvent pool for users to choose from.
 
 Solvents appear on this list will be considered as candidates of multi-solvent systems.

 <p>
  <img src="https://github.com/xueannafang/hsp-toolkits/blob/main/figs/input_solv_sel_example.png">
 </p> 


   Users can edit this form depending on their preference.
   
   *To remove an entry:*
   
   - Select the whole row -> right click -> delete.
   - Update the number in "Group" column by dragging the first two cells to the last filled cell with "filled series" option.
   <p>
 <img src="https://github.com/xueannafang/hsp-toolkits/blob/main/figs/to_rm_from_pool.png">
 </p>
 
   *Note: The "Group" number must be continuous integers starting from 1 to the total number of solvents.*
   
   *Empty cell, discontinuous or overflowed numbers, or wrong sequence in the "Group" column will lead to error.*
   
   <p>
   <img src="https://github.com/xueannafang/hsp-toolkits/blob/main/figs/ip_solv_sel_cm_er.png" width=500>
   </p>
   
   *To add an entry*
   
   - First search the CAS No. of the solvent in the database (db.xlsx).
   - Search this CAS No. in the solvent candidate list (input_solv_sel.xlsx).
   
   **If it is found:**
   
   The solvent you are looking for is already on the list.
   
   **Do not add it again.** (Repeated entries will cause crash.)
   
   **If it is not found:**
   
   Copy the first three cells (No., CAS, Name) in the database (db) spreadsheet to the corresponding cells (No_db, CAS, Solvent) in the candidate list.
   
   (left: db.xlsx) -> (right: input_solv_sel.xlsx)

   No. -> No_db

   CAS -> CAS

   Name -> Solvent
   
   Note that the "No_db" column is not a must to be filled, but is still recommended to do so for development purpose. (Removing this whole column will however fail the calculation.)
   
   - Then update the "Group" column.
   
   <p>
   <img src="https://github.com/xueannafang/hsp-toolkits/blob/main/figs/to_add_to_pool.png" width=1000>
   </p>
   
- **Database** (db.xlsx)

This database is adapted from [C. M. Hansen, Hansen solubility parameters (HSP), Second edn, 2011, vol. 118](https://www.taylorfrancis.com/books/mono/10.1201/9781420006834/hansen-solubility-parameters-charles-hansen).

D, P, H stands for dispersion, dipolar and hydrogen bond sub-parameters (usually written as delta_D, delta_P, detla_H). The unit is MPa<sup>1/2</sup>.

Users can add alias/abbreviations for solvents with long names, e.g., acetonitrile -> ACN, in the "alias" column.

*Note that in this beginner-friendly version, content not in the first six columns (No., CAS, Name, D, P, H) will not affect the calclation. Only to make it convenient for users to search the request entries.*

*Experienced users can personalise the additional parameters by operating other columns of this database file and include them in DataFrame settings.*


### Step 3 Specify arguments in **mix_pred()** function:

```
mix_pred(self, n = 2, rep_time = 50, std = 0.1, tol_pred = 1, red_tol = 0.01)
```

**mix_pred()** contains six arguments: *self, n, rep_time, std, tol_pred, red_tol*.

Please do not change the sequence of them or remove any item. (The first one, **self**, does not need any modification - just leave it there.)

**n**
    
- maximum number of solvents in each multi-solvent combination
- default = 2
- positive integer.

Note that by meaning the "maximum number of solvents", it is because *Solvent Predictor* has included a validation -> filtration process to make sure the predicated concentration of each solvent is paractical enough to be used by experimental researchers.

For example, a prediction of 0.01% solvent A and 99.99% Solvent B can be regarded as redundant information and can be directly replaced by neat Solvent B.

That means even if a user set **n** = 2 in this case, *Sovlent Predictor* will still return one solvent in the case above.

This tolerance (we call it "tolerance of redundant solvent concentration") can be either specified by users in the "**red_tol**" argument, or default set as 1%.

**target HSP**

- contains three parameters in the order of delta_D, delta_P, delta_H
- Each of them must be a float.

**rep_time**

- repeated calculation times for perturbation applied on the target HSP matrix
- default = 50
- positive integer

**std**

- standard deviation of the gaussian random variable applied for perturbation
- default = 0.1

**tol_pred**

- tolerance of error
- can be regarded as concentration deviated from the target
- Calculated concentration deviate more than this value will be filtered.

**red_tol**

- tolerance of redundant solvent concentration
- Calculated conecentration below than this value will be filtered.


### Step 4 Upload and run

Upload solvent candidate list (input_solv_sel.xlsx) and database (db.xlsx) into **SolvPred()** and execute the calculation cell:

```
class SolvPred():
    
    def __init__(self, input_solv, db):
        self.pred = HSP.SolvPredictor(input_solv, db)

    def mix_pred(self, n = 2, rep_time = 50, std = 0.1, tol_pred = 1, red_tol = 0.01):
        self.pred.run_all(n, 18.0, 1.4, 2.0, rep_time = rep_time, std = std, tol = tol_pred, red_tol = red_tol)

sp = SolvPred(r'input_solv_sel.xlsx', r'db.xlsx')
sp.mix_pred()
```

Please note, if you did not change the name of the two input spreadsheets, you can directly run the above code.

Otherwise you need to update the file name in

```
SolvPred(r'input_solv_sel.xlsx', r'db.xlsx')
```

Here the first arugement corresponds to the solvent candidate list, and the second one is the database.

## Output examples

The calculation normally will be done within 1 minute. Once it has been finished, a folder named "solv_sel" will be created under current working directory (i.e., the place where this toolkit is run). Four output files, including two checkpoint excel spreadsheets, one log file and one final result spreadsheet, will be saved in the corresponding folder.

<p>
 <img src=https://github.com/xueannafang/hsp-toolkits/blob/main/figs/op_dir_solvp.png>
 </p>

For general users, the final result (\*final_result.xlsx) and log of input parameters (\*log_SolvPredict.txt) would be of interest.

The following example is an output of dual-solvent predcition (n=2) with target HSP equals to (18.0, 1.4, 2.0)

- Note that this is the HSP of toluene, which also suggests one of the potential important applications of *Solvent Predictor*, that is to find solvent replacement. 

Let's look at those output files:

### solv_sel_log_SolvPredict.txt

In the log file, the input parameters as discussed in the previous section are stored.

<p>
 <img src=https://github.com/xueannafang/hsp-toolkits/blob/main/figs/op_log_toluene_rep.png>
 </p>

The first line: "Solvent amount", equals to the input **n**.

The fourth to seventh line are the original target HSP, which in this case is D = 18.0, P = 1.4, H = 2.0

The next section is set up of statistical perturbation applied on the target HSP.

In this case, we used the default setting (t, std) = (50, 0.1), meaning that the calculation was carried out for 50 times, with each time a different Gaussian random variable with standard deviation equals to 0.1 applied to each element of the original target HSP.

The final part is the set up for results filtration, where we specified the tolerance of error (**tol_pred**) = 1, meaning any HSP deviated more than 1 from the target will be filtered out; and the tolerance of redundant solvents (**red_tol**) remained as 0.01, meaning that concentration below 1% will be regarded as invalid.

Note that the *tol_pred* and *red_tol* may need adjustment if your output file contains no results.



## Principles

See ref X (Manuscript in progress)


