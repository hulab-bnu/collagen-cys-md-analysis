# Complete NC local interaction panel
# Load NC.pse, then run this script.

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

python
for obj in cmd.get_names('objects'):
    if cmd.get_type(obj) == 'object:molecule':
        cmd.enable(obj)
python end

hide everything, all
show cartoon, all
set cartoon_transparency, 0.00, all
set cartoon_loop_radius, 0.25
set cartoon_sampling, 14

set_color NC_A, [0.55, 0.05, 0.12]
set_color NC_B, [0.76, 0.16, 0.22]
set_color NC_C, [0.94, 0.43, 0.43]
set_color NC_D, [0.00, 0.29, 0.36]
set_color NC_E, [0.00, 0.49, 0.52]
set_color NC_F, [0.28, 0.72, 0.68]
color NC_A, chain A
color NC_B, chain B
color NC_C, chain C
color NC_D, chain D
color NC_E, chain E
color NC_F, chain F

# NC contact: Cys70 in chain C and Cys32 in chain F.
select local_cys, (chain C and resi 70) or (chain F and resi 32)
select local_shell, byres (local_cys expand 8.0)
# Surface uses a wider residue-complete envelope than the cartoon context.
# This prevents the molecular surface from ending abruptly at the panel edge.
select local_surface, byres (local_cys expand 11.0)
set cartoon_transparency, 0.70, not local_shell
set cartoon_transparency, 0.10, local_shell

set_color NC_surface_warm, [0.96, 0.79, 0.79]
set_color NC_surface_cool, [0.72, 0.90, 0.90]
hide surface, all
show surface, local_surface
set surface_color, NC_surface_warm, chain A+B+C
set surface_color, NC_surface_cool, chain D+E+F
set transparency, 0.68, local_surface
set surface_quality, 2
set ambient, 0.42
set direct, 0.16
set specular, 0.08
set shininess, 10

hide sticks, all
hide spheres, all
set_color cys_c_green, [0.10, 0.55, 0.20]
set_color cys_n_blue, [0.10, 0.18, 0.62]
set_color cys_o_red, [0.90, 0.05, 0.08]
set_color cys_s_gold, [0.94, 0.65, 0.08]
show spheres, resn CYS
set sphere_scale, 0.32, resn CYS
set sphere_transparency, 0.00, resn CYS
set sphere_color, cys_c_green, resn CYS and elem C
set sphere_color, cys_n_blue, resn CYS and elem N
set sphere_color, cys_o_red, resn CYS and elem O
set sphere_color, cys_s_gold, resn CYS and elem S
rebuild resn CYS, spheres

show sticks, resn CYS
set stick_radius, 0.13, resn CYS
set stick_color, cys_c_green, resn CYS and elem C
set stick_color, cys_n_blue, resn CYS and elem N
set stick_color, cys_o_red, resn CYS and elem O
set stick_color, cys_s_gold, resn CYS and elem S

hide dashes, all
delete NC_SS_contact
distance NC_SS_contact, (chain C and resi 70 and name SG), (chain F and resi 32 and name SG)
set dash_color, cys_s_gold, NC_SS_contact
set dash_width, 3.5, NC_SS_contact
set dash_length, 0.12, NC_SS_contact
set dash_gap, 0.30, NC_SS_contact

label local_cys and name CA, ""
hide labels, local_cys
set label_color, black
set label_outline_color, -1
set label_size, 30
set label_font_id, 7
# Applies only to the distance label because all residue labels are hidden.
# Shift it away from the Cys spheres for a cleaner local-panel annotation.
set label_position, [3.0, 2.0, 0.0]

hide everything, hydro
zoom local_surface, 7

# png NC_local_complete.png, 1800, 1600, dpi=300, ray=1
