[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/EnzymeML/Giess_2026/HEAD)

# This repository accompanies the paper Giess 2026

> "FAIR NMR: Acquisition, Processing, and Modeling of Time-Resolved NMR Data According to the FAIR Principles"  
> *Analytical Chemistry* (2026), DOI: [10.1021/acs.analchem.6c02376](https://doi.org/10.1021/acs.analchem.6c02376)

The folder structure follows the structure of the paper.  
Here you can find all data and Jupyter Notebooks.

## Contents

**Scenario1**: <sup>1</sup>H-NMR kinetic analysis of a Pd(0)-catalyzed cycloisomerization

- `cycloisomerization_analysis.ipynb`: processing of the raw NMR data with NMRPy
- `cycloisomerization_modeling.ipynb`: kinetic modeling of the resulting time courses
- `modules/`: helper functions used by the analysis notebook
- `data/input/`: raw Bruker FIDs and the EnzymeML template
- `data/output/`: the shipped result (EnzymeML documents, fitted model, concentration tables as CSV, NMRPy data models)
- `data/state/`: saved NMRPy states of the three measurements
- `data/comparison/`: the same time courses integrated with MestReNova, and their fit
- `data/run/`: created when the notebooks are run; holds their results and is not part of the repository

**Scenario2**: <sup>31</sup>P-NMR kinetic analysis of an enzymatic cascade reaction

- `cascade_analysis.ipynb`: processing of the raw NMR data with NMRPy
- `cascade_fit.ipynb`: kinetic modeling of the resulting time courses
- `data/input/`: raw Varian FIDs and the EnzymeML document
- `data/output/`: the EnzymeML document with the concentration time courses

## Citation

If you use the data or notebooks of this repository, please cite the paper:

Gieß, T. R.; Matenchi, K.; Wimmer, E. J.; Range, J. P.; Estes, D.; Rohwer, J.; Pleiss, J. FAIR NMR: Acquisition, Processing, and Modeling of Time-Resolved NMR Data According to the FAIR Principles. *Anal. Chem.* **2026**. https://doi.org/10.1021/acs.analchem.6c02376

## Using the Jupyter Notebooks

### With Binder

You can use Binder to run the Jupyter Notebooks online, by clicking on this badge:  
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/EnzymeML/Giess_2026/HEAD)  
(NOTE: you will be directed to Binder, this will not open a new tab! You can use right-click -> Open in new Tab)  
  
This will start a JupyterLab session in your browser.  
On the left you find the folder structure. Double-click on the folder you want to open, e.g. Scenario2. The Jupyter Notebooks have a small orange image icon. Double-click on the notebook you want to run.  
Now you can interact with the notebook. To run all cells click the icon with the two arrows, left of Download, and click on Restart in the popup.  
Each cell will automatically run from top to bottom. When a code cell was run a number appears in the brackets on the left.  

### Locally

To run the Jupyter Notebooks locally, you need to have Python 3.12, 3.13 or 3.14 installed on your computer.

You can install the required packages using pip. Open a terminal and navigate to the folder where you have cloned this repository, then run:

```bash
pip install -r requirements.txt
```

After installing the required packages, you can start JupyterLab by running:

```bash
jupyter lab
```

This will open JupyterLab in your default web browser. You can then navigate to the folder structure on the left and open the Jupyter Notebooks as described in the Binder section.
