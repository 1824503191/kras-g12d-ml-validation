# KRAS G12D ML Screening Pipeline Validation

This repository contains code and processed analysis outputs used for the MSc dissertation project:

**Validation of an ML Screening Pipeline for KRAS G12D Ligands**

The project evaluates whether an existing Random Forest barrier-prediction model for KRAS G12D ligand screening can transfer reliably to improved ligand-in-QM/MM reference barriers.

The repository documents the analysis workflow used in the dissertation, including representative ligand subset analysis, structural visualisation, external validation, feature importance and variance diagnosis, and reduced-feature high/low transfer testing.

---

## Repository Structure

```text
analysis_batch0092/
structure_analysis/
qmmm_validation_47_complete/
high_low_transfer_reduced_features/
README.md
requirements.txt
.gitignore
```

---

## Folder Overview

### `analysis_batch0092/`

Representative ligand subset analysis from the Myriad screening outputs.

This folder contains scripts and processed outputs for:

- reconstructing matched ligand-feature tables;
- defining low/high RF-predicted barrier groups;
- comparing protein–protein, protein–ligand and protein–GTP distance features;
- running exploratory low/high classifiers;
- calculating feature effect sizes;
- identifying important protein–protein distance features;
- generating figures and PyMOL visualisation inputs.

Main outputs are stored in:

```text
analysis_batch0092/figures/
analysis_batch0092/tables/
analysis_batch0092/pymol_scripts/
analysis_batch0092/pymol_images/
```

---

### `structure_analysis/`

Structural analysis and PyMOL visualisation files.

This folder contains:

- scripts for analysing barrier values;
- scripts for selecting representative low/high barrier ligands;
- PyMOL scripts;
- selected ligand lists;
- processed barrier tables;
- final visualisation figures.

Large molecular structure folders and PyMOL session files are not included.

---

### `qmmm_validation_47_complete/`

Validation of the old Random Forest model against improved ligand-in-QM/MM barriers.

This folder contains:

- validation barrier data;
- prediction outputs;
- validation figures;
- scripts for predicted-versus-computed barrier comparison;
- scripts for feature importance and variance diagnosis.

Large trained model files are not included.

---

### `high_low_transfer_reduced_features/`

Reduced-feature high/low transfer test.

This folder contains:

- scripts for constructing the reduced feature set;
- scripts for training and evaluating high/low classifiers;
- result tables;
- figures used in the transfer analysis.

Large model files, logs and raw training data are not included.

---

## Dependencies

The scripts were mainly written in Python 3.

Main Python packages:

```text
numpy
pandas
scipy
scikit-learn
matplotlib
```

Some structural visualisation files require PyMOL.

To install the main Python dependencies:

```bash
pip install -r requirements.txt
```

---

## Reproducibility Notes

The original analysis was performed on UCL computing resources, mainly Myriad and Bohr. Some scripts may contain absolute paths or assumptions specific to those environments.

The full raw molecular structure datasets, QM/MM outputs, trained Random Forest model files, full feature tables and large intermediate files are not included due to file size and project data restrictions. Therefore, the full computational workflow cannot be reproduced from scratch using this repository alone.

This repository is intended to document:

- the analysis code used in the dissertation;
- the main processed outputs;
- the figures and tables used for reporting;
- the computational order of the analysis workflow.

Users with access to the original project data and excluded model files can rerun the scripts after adjusting paths to their own environment.

---

## Suggested Reproduction Workflow

The analysis can be followed in five main stages.

---

## Stage 1: Representative Ligand Subset Analysis

Folder:

```bash
analysis_batch0092/
```

Purpose:

This stage analyses a representative ligand subset from the Myriad screening outputs. It reconstructs matched ligand-feature tables, defines low/high RF-predicted barrier groups, compares feature groups and generates outputs for structural interpretation.

Main inputs:

```text
Legacy RF-predicted barrier outputs
Model-used distance features
Ligand identifiers
Representative ligand subset files
```

Main scripts:

```bash
python build_batch0092_matched_table.py
python build_batch0092_matched_table_all_jobs.py
python check_model_feature_types.py
python classify_low_high_feature_sets.py
python classify_full_batch_with_baseline_and_selection.py
python compare_batch0092_low_high_features.py
python analyse_top_PROPRO_effects.py
python find_selected_PROPRO_features.py
python map_batch0092_features_to_residues.py
python find_best_pdb_for_mapping.py
python add_cohens_d_ci_for_pro_lig.py
```

Main outputs:

```text
analysis_batch0092/tables/
analysis_batch0092/figures/
analysis_batch0092/pymol_scripts/
analysis_batch0092/pymol_images/
```

Expected results:

- matched feature tables;
- feature group classification summaries;
- low/high barrier feature comparisons;
- selected protein–protein distance features;
- figures used for structural interpretation.

Note:

Some full feature tables are excluded from the repository because of file size restrictions. Therefore, scripts that require the full original feature table may need access to the original project data.

---

## Stage 2: Structural Visualisation

Folder:

```bash
structure_analysis/
```

Purpose:

This stage selects representative low/high barrier ligands and generates PyMOL visualisation files for structural interpretation.

Main inputs:

```text
Barrier summary files
Selected ligand lists
Available ligand structure files
KRAS receptor structure
```

Main scripts:

```bash
python analyse_barriers.py
python collect_selected_pdbqt.py
python collect_structures_from_shared_source.py
python make_pymol_script.py
python make_pymol_low_high.py
```

PyMOL visualisation:

```bash
pymol pymol_low_high.pml
```

Main outputs:

```text
selected_ligands.tsv
low_barrier_ligands.txt
high_barrier_ligands.txt
barrier_values_cleaned.tsv
barrier_bimodal_distribution.png
low_high_ligands_overview.png
pymol_low_high.pml
```

Note:

Large structure folders such as receptor files and selected ligand structure folders are not included. The repository includes scripts, selected ligand lists and final visualisation outputs where possible.

---

## Stage 3: External Validation against Improved Ligand-in-QM/MM Barriers

Folder:

```bash
qmmm_validation_47_complete/
```

Purpose:

This stage validates the old Random Forest model against improved ligand-in-QM/MM reference barriers. It compares predicted barriers with computed barriers and reports correlation, error metrics and prediction compression.

Main inputs:

```text
barriers.csv
validation ligand IDs
model-required features
old Random Forest model file
```

Typical script order:

```bash
cd qmmm_validation_47_complete
python scripts/01_check_complete_47_inputs.py
python scripts/02_predict_old_rf_on_complete_current.py
python scripts/03_stats_plots_complete_current.py
```

Main outputs:

```text
qmmm_validation_47_complete/predictions/
qmmm_validation_47_complete/figures/
```

Expected results:

- predicted barrier values;
- predicted-versus-computed barrier comparison;
- Pearson and Spearman correlations;
- confidence intervals;
- MAE and RMSE;
- validation scatter plot.

Note:

The trained Random Forest `.pkl` model file is excluded from this repository due to file size and project data restrictions. Therefore, the prediction script requires access to the original model file before it can be rerun.

---

## Stage 4: Feature Importance and Variance Diagnosis

Folder:

```bash
qmmm_validation_47_complete/
```

Purpose:

This stage compares Random Forest feature importance with feature variance across the improved ligand-in-QM/MM validation structures. It diagnoses whether prediction compression is caused by important features being fixed or near-zero-variance.

Main inputs:

```text
model-used feature list
feature values across validation structures
Random Forest feature importance
validation structures or processed feature outputs
```

Typical script:

```bash
cd qmmm_validation_47_complete
python scripts/04_feature_importance_variance_complete_current.py
```

Main outputs:

```text
qmmm_validation_47_complete/figures/
qmmm_validation_47_complete/model_features/
```

Expected results:

- feature importance versus variance analysis;
- summary of importance on fixed or near-zero-variance features;
- feature group importance for PRO-PRO, PRO-GTP and PRO-LIG features;
- diagnostic figure used in the dissertation.

Note:

Large trained model files and some full intermediate feature files are excluded, so rerunning this stage may require access to original project data.

---

## Stage 5: Reduced-Feature High/Low Transfer Test

Folder:

```bash
high_low_transfer_reduced_features/
```

Purpose:

This stage tests whether the old low/high barrier signal transfers to improved ligand-in-QM/MM labels after removing fixed or near-zero-variance features.

Main inputs:

```text
old screening feature table
new ligand-in-QM/MM validation labels
reduced feature list
features varying across the validation structures
```

Typical script order:

```bash
cd high_low_transfer_reduced_features
python scripts/01_prepare_reduced_features_and_new_labels.py
python scripts/02_train_old_test_new_transfer.py
```

Main outputs:

```text
high_low_transfer_reduced_features/tables/
high_low_transfer_reduced_features/figures/
```

Expected results:

- reduced feature set composition;
- logistic regression and Random Forest transfer-test results;
- ROC AUC values;
- confidence intervals;
- permutation p-values;
- summary tables used in the dissertation.

Note:

Large raw training data, trained model files and logs are excluded from the repository. The included scripts and processed outputs document the transfer-test workflow.

---

## Excluded Files

The following types of files are excluded from the repository:

- trained Random Forest `.pkl` model files;
- full `.parquet` feature tables;
- full molecular structure folders;
- raw docking or structure files;
- QM/MM output folders;
- server log files;
- PyMOL session files;
- large intermediate files.

Examples of excluded files include:

```text
rf_model_dist_below_10_new.pkl
batch0092_feature_table.parquet
high_pdbqt/
low_pdbqt/
receptor/
logs/
models/
old_data/
pymol_sessions/
pymol_local_package/
```

These files are excluded due to file size and project data restrictions.

---

## Notes for Assessors

This repository is provided to support assessment of the code and analysis workflow behind the dissertation.

Because some large input data and trained model files are excluded, not every script is expected to run end-to-end from this repository alone. The repository instead provides the code structure, processed outputs, figures and analysis workflow used to generate the dissertation results.

---
