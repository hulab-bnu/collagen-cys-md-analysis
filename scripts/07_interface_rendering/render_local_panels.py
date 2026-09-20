import os
from pathlib import Path
from pymol import cmd

REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE = Path(os.environ.get("COLLAGEN_PYMOL_SESSIONS", REPO_ROOT / "data" / "pymol_sessions"))
OUTPUT = Path(os.environ.get("COLLAGEN_RESULTS_ROOT", REPO_ROOT / "results")) / "local_panels_png"
SCRIPTS = REPO_ROOT / "pymol"

JOBS = [
    (SOURCE / "NC.pse", SCRIPTS / "NC_local_complete.pml", OUTPUT / "NC_C70_F32_local.png"),
    (SOURCE / "CC.pse", SCRIPTS / "CC_local_C53_Dminus1_Eminus1.pml", OUTPUT / "CC_C53_Dminus1_Eminus1_local.png"),
    (SOURCE / "CC.pse", SCRIPTS / "CC_local_A20_F32.pml", OUTPUT / "CC_A20_F32_local.png"),
]

OUTPUT.mkdir(exist_ok=True)

for session, style, output in JOBS:
    cmd.reinitialize()
    cmd.load(str(session), quiet=1)
    cmd.do("@ " + str(style).replace("\\", "/"))
    cmd.set("ray_opaque_background", 0)
    draw_output = output.with_name(output.stem + "_draw.png")
    cmd.png(str(draw_output), 2400, 1800, dpi=300, ray=0, quiet=1)
    print(draw_output)

cmd.quit()
