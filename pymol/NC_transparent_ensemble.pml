reinitialize
load "data/ensemble_dccm/structures/NC_ensemble.pdb", ensemble
load "data/ensemble_dccm/structures/NC_putty.pdb", reference
set all_states, on
set ray_opaque_background, off
set ray_shadows, off
set antialias, 2
set orthoscopic, on
set cartoon_sampling, 14
set_color ensemble_grey, [0.69, 0.78, 0.80]
set_color chain_a, [0.18, 0.49, 0.66]
set_color chain_b, [0.26, 0.66, 0.60]
set_color chain_c, [0.52, 0.77, 0.58]
set_color cysteine_gold, [0.93, 0.60, 0.12]
hide everything, all
show cartoon, ensemble
color ensemble_grey, ensemble
set cartoon_transparency, 0.78, ensemble
show cartoon, reference
cartoon putty, reference
set cartoon_putty_transform, 0, reference
set cartoon_putty_scale_min, 0.55, reference
set cartoon_putty_scale_max, 1.25, reference
color chain_a, reference and chain A
color chain_b, reference and chain B
color chain_c, reference and chain C
show spheres, reference and resn CYS+CYX+CYM
show sticks, reference and resn CYS+CYX+CYM
color cysteine_gold, reference and resn CYS+CYX+CYM
set sphere_scale, 0.25, reference and resn CYS+CYX+CYM
set stick_radius, 0.10, reference and resn CYS+CYX+CYM
orient reference
zoom reference, 1.18
png results/ensemble_dccm/NC_transparent_ensemble.png, 2200, 560, 300, 1
quit
