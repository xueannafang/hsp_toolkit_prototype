Last update: 01/11/2022

# HSP-toolkits (v1.0)

* *Solvent Predictor*: Based on the target Hansen solubility parameters (HSPs), propose a list of multi-solvent combination.
* *M Locator*: Predict HSPs of the studied material based on a solubility score.

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

*Python 3.7, Juypter Notebook, Microsoft Excel*

The following python packages are required to be installed:

*numpy, pandas, scipy, itertools, abc, os*

## *Solvent Predictor*

Given that HSP of a solvent mixture follows a linear combination of each individual component, researchers can easily calculate the HSP of solvent systems with known components. It is however much more difficult to reverse this process.

The aim of *Solvent Predictor* is to support the solvent suggestion when a desired goal of HSP is known.

The key function is to convert the target HSP into a multi-solvent list based on the requirement of user.

# Set up

* Download the folder of **HSP_SolventPredictor** to local working directory on Windows.
* Open **Solv_pred_class.ipynb** using Juypter Notebook.
* ** Note: **HSP_SolvP.py** contains key calculation process for *Solvent Predictor*. It is imported at the first step. Please be careful to change its name.**
 
# Run the code

Step 1: Execute the first cell (import all the related packages). 

```
import numpy as np
import pandas as pd
from scipy.linalg import pinv
import itertools
import abc
import HSP_SolvP as HSP
import os
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
