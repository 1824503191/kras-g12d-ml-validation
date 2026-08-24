# KRAS G12D ML Screening Pipeline Validation

This repository contains code and processed analysis outputs used for the MSc dissertation project:

**Validation of an ML Screening Pipeline for KRAS G12D Ligands**

The project evaluates an existing machine-learning-assisted screening pipeline for KRAS G12D ligands. The main aim is to test whether an old Random Forest barrier-prediction model can reliably transfer to improved ligand-in-QM/MM reference barriers.

## Repository Structure

### `analysis_batch0092/`

This folder contains scripts and processed outputs for the representative ligand subset analysis from the Myriad screening outputs.

It includes code for:

- reconstructing matched ligand-feature tables;
- comparing low- and high-barrier groups;
- analysing protein–protein, protein–ligand and protein–GTP distance features;
- calculating effect sizes;
- generating figures and PyMOL visualisation inputs.

This folder is related mainly to the representative ligand subset analysis and structural interpretation in the dissertation.

### `structure_analysis/`

This folder contains structural analysis and PyMOL visualisation files.

It includes:

- scripts for selecting low- and high-barrier ligands;
- PyMOL script files;
- selected ligand lists;
- processed barrier tables;
- final structural visualisation figures.

Large molecular structure folders and PyMOL session files are not included.

### `qmmm_validation_47_complete/`

This folder contains scripts and processed outputs for validating the old Random Forest model against improved ligand-in-QM/MM reference barriers.

It includes:

- validation ligand information;
- prediction outputs;
- validation figures;
- scripts for comparing predicted and computed barriers;
- feature importance and variance diagnosis scripts.

Large trained model files and full intermediate feature files are not included.

### `high_low_transfer_reduced_features/`

This folder contains scripts, tables and figures for the reduced-feature high/low transfer test.

It includes:

- scripts for constructing the reduced feature set;
- transfer-test scripts;
- result tables;
- figures used to support the dissertation analysis.

Large model files, logs and raw training data are not included.

## Data Availability

The full molecular structure datasets, QM/MM outputs, trained model files and large intermediate files are not included due to file size and project data restrictions.

In particular, large files such as trained Random Forest `.pkl` models and full feature tables are excluded from this repository. The repository is intended to document the analysis workflow, scripts, processed outputs and figures used in the dissertation.

## Main Analysis Components

The repository supports the following parts of the dissertation:

1. Large-scale screening output analysis.
2. Representative ligand subset feature analysis.
3. Structural interpretation and PyMOL visualisation.
4. Validation of the old Random Forest model against improved ligand-in-QM/MM barriers.
5. Feature importance and variance diagnosis.
6. Reduced-feature high/low transfer testing.

## Author

Mingsong Geng  
MSc Scientific and Data Intensive Computing  
University College London
