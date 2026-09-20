# ================================================================
# Complete mC local interaction panel
#
# Step 1: Load mC.pse in PyMOL.
# Step 2: Run this file with File > Run Script.
# This script expects chains A to F.
# ================================================================

# ---------- 1. Stable white background and transparency ----------
bg_color white
set orthoscopic, on
set depth_cue, off
set antialias, 2
set ray_shadows, off
set ray_opaque_background, off
set ray_interior_color, white
set ray_interior_shadows, off
set use_shaders, on
set transparency_mode, 3
set backface_cull, off
set two_sided_lighting, -1

# Ensure the molecular object stored in the PSE session is enabled.
python
for obj in cmd.get_names('objects'):
    if cmd.get_type(obj) == 'object:molecule':
        cmd.enable(obj)
python end

# ---------- 2. Cartoon background ----------
hide everything, all
show cartoon, all
set cartoon_transparency, 0.00, all
set cartoon_loop_radius, 0.25
set cartoon_sampling, 14

# Wine red and peacock teal chain palette.
set_color mC_A, [0.55, 0.05, 0.12]
set_color mC_B, [0.76, 0.16, 0.22]
set_color mC_C, [0.94, 0.43, 0.43]
set_color mC_D, [0.00, 0.29, 0.36]
set_color mC_E, [0.00, 0.49, 0.52]
set_color mC_F, [0.28, 0.72, 0.68]
color mC_A, chain A
color mC_B, chain B
color mC_C, chain C
color mC_D, chain D
color mC_E, chain E
color mC_F, chain F

# ---------- 3. Contact region and its surface ----------
# mC contact pair: Cys37 in chain A and Cys16 in chain E.
select local_cys, (chain A and resi 37) or (chain E and resi 16)
select local_shell, byres (local_cys expand 8.0)
# Use a wider, residue-complete envelope for the surface so its edge is smooth.
select local_surface, byres (local_cys expand 11.0)

# Fade remote cartoon and keep the local backbone stronger.
set cartoon_transparency, 0.70, not local_shell
set cartoon_transparency, 0.10, local_shell

# The transparent surface covers all amino acids in the local environment.
# Its colors are independent from Cys spheres and sticks.
set_color mC_surface_warm, [0.96, 0.79, 0.79]
set_color mC_surface_cool, [0.72, 0.90, 0.90]
hide surface, all
show surface, local_surface
set surface_color, mC_surface_warm, chain A+B+C
set surface_color, mC_surface_cool, chain D+E+F
set transparency, 0.68, local_surface
set surface_quality, 2
set ambient, 0.42
set direct, 0.16
set specular, 0.08
set shininess, 10

# ---------- 4. All Cys residues as CPK ball and stick ----------
# All Cys atoms use a common sphere size of 0.32.
# Carbon is green, nitrogen blue, oxygen red, sulfur gold.
hide sticks, all
hide spheres, all
set_color cys_c_green, [0.10, 0.55, 0.20]
set_color cys_n_blue, [0.10, 0.18, 0.62]
set_color cys_o_red, [0.90, 0.05, 0.08]
set_color cys_s_gold, [0.94, 0.65, 0.08]
show spheres, resn CYS
set sphere_scale, 0.32, resn CYS
set sphere_scale, 0.32, resn CYS and elem S
set sphere_transparency, 0.00, resn CYS
set sphere_color, cys_c_green, resn CYS and elem C
set sphere_color, cys_n_blue, resn CYS and elem N
set sphere_color, cys_o_red, resn CYS and elem O
set sphere_color, cys_s_gold, resn CYS and elem S
rebuild resn CYS, spheres

# Overlay thin CPK sticks on the Cys spheres.
show sticks, resn CYS
set stick_radius, 0.13, resn CYS
set stick_color, cys_c_green, resn CYS and elem C
set stick_color, cys_n_blue, resn CYS and elem N
set stick_color, cys_o_red, resn CYS and elem O
set stick_color, cys_s_gold, resn CYS and elem S

# ---------- 5. S to S contact and distance label ----------
hide dashes, all
delete mC_SS_contact
distance mC_SS_contact, (chain A and resi 37 and name SG), (chain E and resi 16 and name SG)
set dash_color, cys_s_gold, mC_SS_contact
set dash_width, 3.5, mC_SS_contact
set dash_length, 0.12, mC_SS_contact
set dash_gap, 0.30, mC_SS_contact

# Keep only the distance number. Do not show Cys residue names.
label local_cys and name CA, ""
hide labels, local_cys
set label_color, black
set label_outline_color, -1
set label_size, 30
set label_font_id, 7
# Move the distance number away from the Cys spheres.
set label_position, [3.0, 2.0, 0.0]

# ---------- 6. Remove hydrogens and frame the local panel ----------
hide everything, hydro
zoom local_surface, 7

# ---------- 7. Transparent PNG export ----------
# png mC_local_complete.png, 1800, 1600, dpi=300, ray=1
