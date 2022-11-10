Last update: 07/11/2022

The aim of this project is to develop quantitative solvent selection tools for synthesis and property control of functional materials. The focus is based on Hansen solubility parameters (HSP).

HSP decompose the molecular interaction into dispersion (𝜹𝑫), dipolar (𝜹𝑷) and hydrogen-bonding (𝜹𝑯) interactions.

Solubility can be quantified by the Hansen distance (R). HSP of a solvent mixture follows a linear combination of its components, which greatly expands options based on given lab solvents.

<p>
 <img src=https://github.com/xueannafang/hsp-toolkits/blob/main/figs/HSP_general.png width=300>
 </p>

There are nevertheless a couple of challenges when using HSP:

1. Determining which solvent combinations can lead to a given HSP.
2. Obtaining HSP of the studied materials (M).

To improve these issues, we developed two python-based toolkits:

* [*Solvent Predictor*](https://github.com/xueannafang/hsp-toolkits/blob/main/HSP_SolventPredictor/solv_pred_readme.md): Based on the target Hansen solubility parameters (HSPs), propose a list of multi-solvent combination.

* [*M Locator*](https://github.com/xueannafang/hsp-toolkits/edit/main/HSP_MLocator/mloc_readme.md): Predict HSPs of the studied material based on an experimental-measurable solubility score.

 <p>
  <img src="https://github.com/xueannafang/hsp-toolkits/blob/main/figs/sch_sp_mloc.png" width=700>
 </p>

These toolkits are developed with Python 3.7.3 and tested in Windows 10.

<p>
 <img src=https://github.com/xueannafang/hsp-toolkits/blob/main/figs/logo_all.png width=200>
 </p>
 
 
 ### References

1. C. Hansen, Hansen Solubility Parameters – A user’s handbook, 2nd edition, 2011.
2. X. Fang, C. F.J. Faul, N. Fey, E. Gale, Solvent Predictor - A python toolkit to predict multi-solvent combinations with target Hansen solubility parameters (manuscript in preparation)
3. X. Fang, U. Karatayeva, S. A. A. Siyabi, B. B. Narzary, M. G. Girgin, D. Mukhanov, C. F.J. Faul, N. Fey, E. Gale, M Locator - A python toolkit to predict Hansen solubility parameters of functional materials (manuscript in preparation).


---

This project is developed by

- [Xue Fang](https://www.linkedin.com/in/xue-fang-811204163/) (School of chemistry, University of Bristol, UK)

with instruction from

- [Prof Charl FJ Faul](https://faulresearchgroup.com/charl-f-j-faul/) (School of chemistry, University of Bristol, UK)
- [Prof Natalie Fey](https://feygroupchem.wordpress.com/) (School of chemistry, University of Bristol, UK)
- [Dr Ella Gale](https://www.bristol.ac.uk/people/person/Ella-Gale-58ab10ba-8b85-4513-944e-6d9020b6ff2c/) (School of chemistry, University of Bristol, UK)

Experimental test data of *M Locator* is contributed by the Conjugated Microporous Polymers (CMP) team of the [Faul research group](https://faulresearchgroup.com/).

The author acknowledges the following organisations for funding this project:

- [University of Bristol](https://www.bristol.ac.uk/)
- [Chinese Scholarship Council](https://www.chinesescholarshipcouncil.com/)
- [Royal Society of Chemistry Researcher Development Grants](https://www.rsc.org/prizes-funding/funding/find-funding/researcher-development-grant/)

The author also thanks all who have provided inputs and technical support during the design and exploration of these toolkits, in particular:

- [Dr Jie Chen](https://scholar.google.com/citations?user=GPM9kTgAAAAJ&hl=en) (Fuzhou University, China)
- [Bo Gao](https://www.linkedin.com/in/bo-gao-771841199/) (School of physics, University of Bristol, UK)


---

This project is licensed under [GPL-3.0](https://www.gnu.org/licenses/gpl-3.0.html).
