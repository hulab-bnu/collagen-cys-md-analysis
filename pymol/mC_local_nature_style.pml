# ================================================================
# mC local interaction panel in a Nature-like structural style
#
# Run after the full overview styling script.
# Focus: Cys37 in chain A and Cys16 in chain E.
# ================================================================

# Store the overview before making the local panel.
scene overview, store

# ---------- 1. Define the contact and its local environment ----------
select local_cys, (chain A and resi 37) or (chain E and resi 16)
select local_shell, byres (local_cys expand 8.0)

# ---------- 2. Keep a restrained structural context ----------
# The local panel uses a transparent surface around the contact region.
# It inherits the overview surface colors while keeping the Cys readable.
hide surface, all
show surface, local_shell
set transparency, 0.68, local_shell
set surface_quality, 1

# The local backbone remains clear while the distant backbone is faded.
show cartoon, all
set cartoon_transparency, 0.70, not local_shell
set cartoon_transparency, 0.10, local_shell

# ---------- 3. Cys contact ----------
# All non Cys residues remain as cartoons in this local panel.
hide spheres, all
hide sticks, all

# All Cys residues are spheres with an optimized CPK palette.
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
hide everything, hydro

# ---------- 4. Contact line and concise labels ----------
hide dashes, all
delete mC_SS_contact
distance mC_SS_contact, (chain A and resi 37 and name SG), (chain E and resi 16 and name SG)
set dash_color, cys_s_gold, mC_SS_contact
set dash_width, 3.5, mC_SS_contact
set dash_length, 0.25, mC_SS_contact
set dash_gap, 0.12, mC_SS_contact

label local_cys and name CA, ""
hide labels, local_cys
set label_color, black
set label_outline_color, white
set label_size, 22

# ---------- 5. Camera and transparent export ----------
zoom local_shell, 7
set ray_opaque_background, off

# png mC_local_nature_style.png, 1800, 1600, dpi=300, ray=1

# Return to the overview after exporting with:
# scene overview, recall
