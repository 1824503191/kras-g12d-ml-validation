reinitialize

load /myriadfs/home/ucapmge/Scratch/omm_pipeline_screening/manual_charmm_setup/step1_pdbreader.pdb, receptor

hide everything
show cartoon, receptor
color grey80, receptor

# Top positive features:
# high-barrier ligands are farther from these residues
select high_farther_region, chain P and resi 60+61+62+67+68+69+88+1235+1283
show sticks, high_farther_region
color yellow, high_farther_region

# Main negative features:
# high-barrier ligands are closer to these residues
select high_closer_region, chain P and resi 1241+1237+65+63
show sticks, high_closer_region
color magenta, high_closer_region

# Add labels for most important residues
label chain P and resi 60 and name CA, "GLY60"
label chain P and resi 69 and name CA, "ASP69"
label chain P and resi 61 and name CA, "GLN61"
label chain P and resi 1241 and name CA, "ARG1241"
label chain P and resi 88 and name CA, "LYS88"

set label_size, 16
set label_color, black
set stick_radius, 0.18

orient receptor
zoom high_farther_region or high_closer_region, 12

bg_color white
set ray_opaque_background, off

png /myriadfs/home/ucapmge/Scratch/omm_pipeline_screening/analysis_batch0092/pymol_images/batch0092_feature_residue_mapping.png, 
width=1800, height=1400, dpi=300, ray=1

save /myriadfs/home/ucapmge/Scratch/omm_pipeline_screening/analysis_batch0092/pymol_sessions/batch0092_feature_residue_mapping.pse
