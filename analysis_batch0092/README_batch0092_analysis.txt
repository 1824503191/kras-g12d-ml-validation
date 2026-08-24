README: batch0092 full analysis package
Prepared by: Mingsong
Date: 2026-07-27

Project purpose
---------------
This folder contains the updated full-batch analysis for batch0092 using the old RF model's predicted barriers. The aim is to understand which structural distance features are associated with low/high RF-predicted barriers, with particular attention to PRO-LIG, PRO-GTP and PRO-PRO feature groups.

Important note
--------------
The barriers analysed in this folder are RF-predicted barriers from the old ML model, not newly computed QM/MM barriers and not experimental values. Therefore, the results should be interpreted as model-interpretation evidence rather than direct physical proof.

Dataset summary
---------------
The full matched batch0092 table contains:

- Matched ligands: 992
- Low RF-predicted barrier group (< 28 kcal/mol): 285
- Middle group (28–29 kcal/mol): 126
- High RF-predicted barrier group (>= 29 kcal/mol): 581
- Low/high samples used for classifier analysis after excluding the middle band: 866

Feature groups
--------------
The old RF model uses 4,086 distance features:

- PRO-PRO residue-residue features: 3,982
- PRO-GTP residue-GTP features: 83
- PRO-LIG residue-ligand features: 21

Main analyses included
----------------------
1. Full matched table construction
   Script:
   - build_full_batch.log
   - build_full_batch0092_matched_table.py
   - build_batch0092_matched_table_all_jobs.py

   Output:
   - tables/batch0092_features_with_barriers.parquet

2. Feature group classifier comparison
   Script:
   - classify_full_batch_with_baseline_and_selection.py

   Aim:
   Compare whether different feature groups can separate low vs high RF-predicted barrier ligands.

   Main result:
   - PRO-LIG features alone showed weak signal.
   - PRO-GTP features showed modest signal.
   - Selected PRO-PRO residue-residue features showed the strongest low/high separation.
   - Majority baseline was included, and balanced accuracy / ROC AUC were used to avoid misleading raw accuracy from class imbalance.

3. PRO-LIG effect-size analysis
   Script:
   - add_cohens_d_ci_for_pro_lig.py

   Output:
   - tables/batch0092_PRO_LIG_cohens_d_with_CI.tsv
   - figures/batch0092_PRO_LIG_cohens_d_with_CI.png/pdf/svg

   Main result:
   Direct ligand-residue distance effects were generally weak. Only PRO60-LIG and PRO61-LIG showed small Cohen's d effects with confidence intervals not crossing zero. Most PRO-LIG differences should be treated as exploratory.

4. Stable PRO-PRO feature selection
   Script:
   - find_selected_PROPRO_features.py

   Output:
   - tables/batch0092_PROPRO_selected_feature_frequency.tsv

   Main result:
   The most consistently selected PRO-PRO features mainly involve GLU62 and residues around 1275-1280, including:
   - GLU62-SER1279
   - PHE1275-ARG1276
   - GLU62-LEU1280
   - THR87-ASN1278
   - VAL14-GLU62
   - THR124-ASN1278
   - PHE90-ASN1278
   - SER89-ASN1278

5. Direction of top PRO-PRO effects
   Script:
   - analyse_top_PROPRO_effects.py

   Output:
   - tables/batch0092_top_PROPRO_feature_effects.tsv
   - figures/batch0092_top_PROPRO_cohens_d_with_CI.png/pdf/svg

   Main result:
   Several top selected PRO-PRO distances are smaller in the high RF-predicted barrier group. This supports the interpretation that the old RF model's low/high predictions are more strongly associated with protein conformational reorganisation than with direct ligand-residue contacts alone.

6. Structural mapping
   Relevant folders:
   - pymol_scripts/
   - pymol_images/
   - pymol_sessions/
   - pymol_local_package/

   Purpose:
   These files are used to map selected residues and representative distance features onto the protein structure in PyMOL.

   Important caveat:
   Dashed distance lines in PyMOL figures should be described as representative visual guides. Unless otherwise confirmed, they should not be described as the exact atom-pair distances used in feature generation.

Folder guide
------------
tables/
  Main numerical outputs, including matched tables, classifier results, Cohen's d tables, selected feature frequency tables, and PRO-PRO effect direction tables.

figures/
  Statistical plots generated from the analysis, including Cohen's d confidence interval plots.

pymol_images/
  Exported PyMOL images.

pymol_sessions/
  PyMOL session files.

pymol_scripts/
  PyMOL scripts used to generate structural figures.

pymol_local_package/
  Files prepared for local PyMOL visualisation.

*.py scripts
  Analysis scripts used for feature grouping, classifier comparison, effect-size analysis, and residue mapping.

*.log files
  Terminal logs from key analysis runs.

Main interpretation
-------------------
The full batch0092 analysis suggests that the old RF model's low/high predicted barrier labels are associated much more strongly with protein residue-residue distance changes than with direct ligand-residue distances. In particular, the repeated selection of GLU62 and residues around 1275-1280 suggests a possible protein conformational reorganisation signal.

This should be presented as an exploratory model-interpretation result, not as causal proof of the true physical mechanism.

Suggested short summary
-----------------------
In the full batch0092 analysis, PRO-LIG ligand-residue distances showed only weak effects, while selected PRO-PRO residue-residue distances strongly separated low and high RF-predicted barrier groups. The most stable PRO-PRO signals involved GLU62 and residues around 1275-1280, supporting the interpretation that the old RF model is capturing protein conformational reorganisation rather than direct ligand-residue contact alone.
