# KRAS G12D ML Screening Pipeline Validation

This repository contains code and processed analysis outputs used for the MSc dissertation project:

**Validation of an ML Screening Pipeline for KRAS G12D Ligands**

The project evaluates an existing machine-learning-assisted screening pipeline for KRAS G12D ligands. The main aim is to assess whether an older Random Forest barrier-prediction model can reliably transfer to improved ligand-in-QM/MM reference barriers.

The repository documents the analysis workflow used in the dissertation, including large-scale screening analysis, representative ligand subset feature analysis, structural visualisation, external validation, feature importance and variance diagnosis, and reduced-feature high/low transfer testing.

---

## Repository Structure

### `analysis_batch0092/`

This folder contains scripts and processed outputs for the representative ligand subset analysis from the Myriad screening outputs.

It includes code and outputs for:

- reconstructing matched ligand-feature tables;
- comparing low- and high-barrier groups;
- analysing protein–protein, protein–ligand and protein–GTP distance features;
- calculating effect sizes such as Cohen's d;
- selecting important protein–protein distance features;
- generating figures and PyMOL visualisation inputs.

This folder mainly supports the representative ligand subset analysis and structural interpretation sections of the dissertation.

---

### `structure_analysis/`

This folder contains structural analysis and PyMOL visualisation files.

It includes:

- scripts for analysing barrier values;
- scripts for selecting low- and high-barrier ligands;
- PyMOL script files;
- selected ligand lists;
- processed barrier tables;
- final structural visualisation figures.

Large molecular structure folders and PyMOL session files are not included.

---

### `qmmm_validation_47_complete/`

This folder contains scripts and processed outputs for validating the old Random Forest model against improved ligand-in-QM/MM reference barriers.

It includes:

- validation ligand information;
- computed barrier tables;
- prediction outputs;
- validation figures;
- scripts for comparing predicted and computed barriers;
- scripts for feature importance and variance diagnosis.

This folder mainly supports the external validation and feature-space diagnosis sections of the dissertation.

Large trained model files and full intermediate feature files are not included.

---

### `high_low_transfer_reduced_features/`

This folder contains scripts, tables and figures for the reduced-feature high/low transfer test.

It includes:

- scripts for constructing the reduced feature set;
- scripts for training and evaluating high/low classifiers;
- result tables;
- figures used to support the reduced-feature transfer analysis.

Large model files, logs and raw training data are not included.

---

## Main Analysis Components

This repository supports the following parts of the dissertation:

1. Large-scale Myriad screening output analysis.
2. Representative ligand subset feature analysis.
3. Structural interpretation and PyMOL visualisation.
4. Validation of the old Random Forest model against improved ligand-in-QM/MM barriers.
5. Feature importance and variance diagnosis.
6. Reduced-feature high/low transfer testing.

---

## Dependencies

The analysis scripts were mainly written in Python 3.

The main Python packages used include:

- numpy
- pandas
- scipy
- scikit-learn
- matplotlib

Some structural visualisation files require PyMOL.

The original analysis was performed on UCL computing resources, including Myriad and Bohr. Some scripts may contain paths or assumptions specific to those computing environments and may need to be adjusted before being run on a different machine.

---

## How to Use This Repository

The full raw molecular datasets, complete QM/MM outputs and trained model files are not included. Therefore, the entire computational workflow cannot be reproduced from scratch using this repository alone.

However, the repository contains the scripts, processed outputs and figures used to document the analysis workflow in the dissertation.

Suggested reading order:

1. `analysis_batch0092/`  
   Representative ligand subset analysis, feature comparison, low/high grouping and structural feature interpretation.

2. `structure_analysis/`  
   Structural visualisation scripts, PyMOL files, selected ligand lists and visual outputs.

3. `qmmm_validation_47_complete/`  
   Validation of the old Random Forest model against improved ligand-in-QM/MM computed barriers, including prediction outputs and validation figures.

4. `high_low_transfer_reduced_features/`  
   Reduced-feature high/low transfer test using only features that varied across the improved ligand-in-QM/MM structures.

---

## Data Availability and Excluded Files

Large files are not included due to file size and project data restrictions.

Excluded files include:

- trained Random Forest `.pkl` model files;
- full `.parquet` feature tables;
- full molecular structure folders;
- QM/MM output folders;
- raw docking or structure files;
- large intermediate files;
- server log files;
- PyMOL session files.

The repository instead provides scripts, processed summary outputs and figures used in the dissertation analysis.

---

## Notes on Reproducibility

This repository is intended to document the code and processed analysis outputs behind the dissertation results.

Because some full input data and trained model files are excluded, users may not be able to run every script directly without access to the original project data. The scripts are included to show the computational procedures used for:

- generating screening summaries;
- analysing structural distance features;
- validating predicted barriers against computed barriers;
- diagnosing feature importance and feature variance;
- testing high/low transferability using reduced features.

---

## Author

Mingsong Geng  
MSc Scientific and Data Intensive Computing  
University College London  
August 2026
