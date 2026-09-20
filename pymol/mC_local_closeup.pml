# ================================================================
# mC local interaction close up
#
# Run this after the full collagen styling script.
# It focuses on the A37 and E16 cysteine contact in mC.
# ================================================================

# Save the current overview scene for later restoration.
scene overview, store

# ---------- 1. Interaction center and local environment ----------
select local_cys, (chain A and resi 37) or (chain E and resi 16)
select local_shell, byres (local_cys expand 6.0)

# ---------- 2. Context: a small transparent local surface ----------
hide surface, all
show surface, local_shell
set transparency, 0.65, local_shell

# Keep the contact region opaque and fade the rest of the cartoon.
show cartoon, all
set cartoon_transparency, 0.62, not local_shell
set cartoon_transparency, 0.00, local_shell

# ---------- 3. Chemistry: Cys as CPK spheres ----------
# Only the two interacting Cys are shown as spheres in this local panel.
hide spheres, all
show spheres, local_cys
set sphere_scale, 0.50, local_cys
set sphere_scale, 0.60, local_cys and elem S
set sphere_transparency, 0.00, local_cys

# These representation specific colors leave the cartoons unchanged.
set sphere_color, cys_carbon, local_cys and elem C
set sphere_color, cys_nitrogen, local_cys and elem N
set sphere_color, cys_oxygen, local_cys and elem O
set sphere_color, cys_sulfur_focus, local_cys and elem S

# Neighboring side chains provide restrained chemical context.
set_color local_neighbor_gray, [0.40, 0.43, 0.46]
show sticks, local_shell and sidechain and not local_cys
set stick_radius, 0.13, local_shell and sidechain and not local_cys
set stick_color, local_neighbor_gray, local_shell and sidechain and not local_cys

# ---------- 4. The S to S contact ----------
delete mC_S_contact
distance mC_S_contact, (chain A and resi 37 and name SG), (chain E and resi 16 and name SG)
set dash_color, cys_sulfur_focus, mC_S_contact
set dash_width, 3.0, mC_S_contact

# Optional labels for a supplementary figure or a close up panel.
# label local_cys and name CA, "%s%s" % (chain,resi)
# set label_color, black
# set label_size, 18

# ---------- 5. Camera and export ----------
zoom local_shell, 8
set ray_opaque_background, off

# png mC_local_closeup.png, 1800, 1800, dpi=300, ray=1

# To return to the overview after exporting, use:
# scene overview, recall
