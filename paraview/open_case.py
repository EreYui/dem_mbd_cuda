"""ParaView GUI startup script for the converted DEM/MBD time series."""

import os
import sys
from pathlib import Path

from paraview.simple import (  # type: ignore
    ColorBy,
    GetAnimationScene,
    GetActiveViewOrCreate,
    Glyph,
    Hide,
    OpenDataFile,
    Render,
    SaveState,
    Show,
)


def read_case_name(script_dir: Path) -> str:
    config = script_dir / "case.dat"
    if config.is_file():
        for raw_line in config.read_text(encoding="utf-8-sig").splitlines():
            line = raw_line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            if key.strip().upper() == "CASE_NAME" and value.strip():
                return value.strip()
    return "DATA"


def find_project_root() -> Path:
    """Locate the repository even when ParaView executes this as a string.

    ParaView 6.0 may not define ``__file__`` for ``--script``.  The launcher
    therefore passes an environment variable; the other candidates keep the
    script usable from ParaView's Python shell and from pvpython.
    """

    candidates = []
    configured = os.environ.get("DEM_MBD_PARAVIEW_PROJECT_ROOT")
    if configured:
        candidates.append(Path(configured))

    script_name = globals().get("__file__")
    if script_name:
        candidates.append(Path(script_name).resolve().parents[1])

    for argument in sys.argv:
        if argument.lower().endswith("open_case.py"):
            candidates.append(Path(argument).resolve().parents[1])

    current = Path.cwd().resolve()
    candidates.extend((current, *current.parents))
    for candidate in candidates:
        root = candidate.resolve()
        if (root / "CMakeLists.txt").is_file() and (root / "Data").is_dir():
            return root
    raise RuntimeError(
        "Cannot locate the DEM-MBD project root. Set "
        "DEM_MBD_PARAVIEW_PROJECT_ROOT before starting ParaView."
    )


PROJECT_ROOT = find_project_root()
SCRIPT_DIR = PROJECT_ROOT / "paraview"
configured_case = os.environ.get("DEM_MBD_PARAVIEW_CASE_DIR")
if configured_case:
    CASE_DIR = Path(configured_case).resolve()
else:
    case_name = read_case_name(SCRIPT_DIR)
    if not case_name or case_name in {".", ".."} or "/" in case_name or "\\" in case_name:
        raise RuntimeError(
            "CASE_NAME in paraview/case.dat must be one directory name under Data."
        )
    CASE_DIR = (PROJECT_ROOT / "Data" / case_name).resolve()

OUTPUT = CASE_DIR / "OutputFile" / "paraview"
PARTICLES = OUTPUT / "particles.pvd"
BODIES = OUTPUT / "bodies.pvd"

print(f"ParaView case directory: {CASE_DIR}")

if not PARTICLES.is_file() and not BODIES.is_file():
    raise RuntimeError(
        f"No ParaView time series found under {OUTPUT}. Run prepare_case.bat first."
    )

view = GetActiveViewOrCreate("RenderView")
view.Background = [0.12, 0.14, 0.18]

if PARTICLES.is_file():
    particle_points = OpenDataFile(str(PARTICLES))
    particle_display = Show(particle_points, view)
    try:
        # GPU-instanced representation: much lighter than constructing a sphere
        # mesh for every particle with the Glyph filter.
        particle_display.Representation = "3D Glyphs"
        particle_display.GlyphType = "Sphere"
        particle_display.ScaleArray = ["POINTS", "radius"]
        particle_display.ScaleFactor = 2.0
    except Exception:
        # Compatibility fallback for older ParaView releases.
        Hide(particle_points, view)
        particle_glyphs = Glyph(
            registrationName="Particles (sphere glyphs)",
            Input=particle_points,
            GlyphType="Sphere",
        )
        particle_glyphs.ScaleArray = ["POINTS", "radius"]
        particle_glyphs.ScaleFactor = 2.0
        particle_glyphs.GlyphMode = "All Points"
        particle_glyphs.OrientationArray = ["POINTS", "No orientation array"]
        try:
            particle_glyphs.GlyphType.ThetaResolution = 8
            particle_glyphs.GlyphType.PhiResolution = 8
        except AttributeError:
            pass
        particle_display = Show(particle_glyphs, view)
        particle_display.Representation = "Surface"
    try:
        ColorBy(particle_display, ("POINTS", "speed"))
        particle_display.RescaleTransferFunctionToDataRange(True, False)
    except Exception:
        particle_display.DiffuseColor = [0.76, 0.62, 0.34]

if BODIES.is_file():
    bodies = OpenDataFile(str(BODIES))
    body_display = Show(bodies, view)
    body_display.Representation = "Surface"
    try:
        ColorBy(body_display, ("CELLS", "body_id"))
        body_display.RescaleTransferFunctionToDataRange(True, False)
    except Exception:
        body_display.DiffuseColor = [0.75, 0.78, 0.84]

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
view.ResetCamera()
Render()

state_path = OUTPUT / "dem_mbd_visualization.pvsm"
try:
    SaveState(str(state_path))
    print(f"ParaView state saved to: {state_path}")
except Exception as error:
    print(f"Warning: could not save ParaView state: {error}")
