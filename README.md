# Bachelor's Thesis Repository

## Overview

This repository contains the code developed for my Bachelor's Thesis (TFG), named "Characterizing epileptic seizures from subcutaneous electroencephalographic recordings using univariate phase-based signal analysis measures".
This thesis focused on the analysis of phase-based metrics for characterizing dynamical changes in both simulated systems and epilepsy recordings.

The repository is divided into two main folders:

* **Jupyter Notebooks**: scripts and analyses related to synthetic signals, including the AR(1) process and the Rössler system.
* **MATLAB**: scripts used for the analysis of patient subcutaneous EEG (sqEEG) recordings and seizure-related dynamics.

---

## Repository Structure

```text
.
├── Jupyter/
│   ├── AR1_tests/
│   ├── Rösslers_tests/
│   
│
├── MATLAB/
│   ├── EEGvisual_mat.md/
│   ├── anova_test.md/
│   ├── channels.md/
│   └── ...
│
└── README.md
```

### Jupyter Folder

This folder contains the notebooks used to study the behaviour of the phase-based metrics in controlled dynamical systems.

---

### MATLAB Folder

This folder contains the MATLAB scripts used for the analysis of patient recordings obtained from the SUBER study.


## Requirements

### Python Environment

The Jupyter notebooks were developed using Python 3 and require:

* NumPy
* SciPy
* Matplotlib
* Pandas
* Jupyter Notebook

Install dependencies with:

```bash
pip install numpy scipy matplotlib pandas jupyter
```

### MATLAB

The patient analysis was developed in MATLAB.

Required toolboxes may include:

* Signal Processing Toolbox
* Statistics and Machine Learning Toolbox

depending on the specific script being executed.

---

## Author

Raquel Benjumea

Bachelor's Thesis (TFG)

Universitat Pompeu Fabra (UPF)
