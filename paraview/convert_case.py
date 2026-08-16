#!/usr/bin/env python3
"""Convert DEM/MBD text snapshots into two ParaView time series.

The converter intentionally uses only the Python standard library.  It can be
run with a normal Python 3 interpreter; ParaView Python modules are not needed.
"""

from __future__ import annotations

import argparse
import base64
import csv
import math
import re
import struct
import sys
from array import array
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


PARTICLE_PATTERN = re.compile(r"^ph\.(\d+)\.bt$")
BODY_PATTERN = re.compile(r"^At\.(\d+)\.bt$")
CONFIG_FILE_NAME = "case.dat"


def read_case_config(script_dir: Path) -> dict[str, str]:
    """Read the shared KEY=VALUE launcher configuration."""

    config: dict[str, str] = {}
    path = script_dir / CONFIG_FILE_NAME
    if not path.is_file():
        return config
    for line_number, raw_line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            fail(f"Invalid {path} line {line_number}: expected KEY=VALUE")
        key, value = line.split("=", 1)
        config[key.strip().upper()] = value.strip()
    return config


def case_directory_from_name(project_root: Path, case_name: str) -> Path:
    """Resolve one safe directory name below the repository Data folder."""

    if (
        not case_name
        or case_name in {".", ".."}
        or "/" in case_name
        or "\\" in case_name
    ):
        fail("CASE_NAME must be one directory name under Data, for example DATA or DATA3")
    return project_root / "Data" / case_name


@dataclass(frozen=True)
class BodyMesh:
    name: str
    points: tuple[tuple[float, float, float], ...]
    connectivity: tuple[int, ...]
    offsets: tuple[int, ...]
    cell_types: tuple[int, ...]


@dataclass(frozen=True)
class BodyState:
    position: tuple[float, float, float]
    quaternion: tuple[float, float, float, float]


def fail(message: str) -> "None":
    raise RuntimeError(message)


def numeric_payload(typecode: str, values: Iterable[float | int]) -> str:
    data = array(typecode, values)
    if typecode == "i" and data.itemsize != 4:
        fail("This Python build does not provide 32-bit array('i').")
    if sys.byteorder != "little":
        data.byteswap()
    raw = data.tobytes()
    if len(raw) >= 2**32:
        fail("A single VTK array exceeds the UInt32 binary-header limit.")
    return base64.b64encode(struct.pack("<I", len(raw)) + raw).decode("ascii")


def write_data_array(
    stream,
    vtk_type: str,
    name: str | None,
    components: int,
    typecode: str,
    values: Iterable[float | int],
    indent: str,
) -> None:
    attributes = [f'type="{vtk_type}"']
    if name:
        attributes.append(f'Name="{name}"')
    if components != 1:
        attributes.append(f'NumberOfComponents="{components}"')
    attributes.append('format="binary"')
    stream.write(f"{indent}<DataArray {' '.join(attributes)}>")
    stream.write(numeric_payload(typecode, values))
    stream.write("</DataArray>\n")


def read_step_size(world_file: Path) -> float:
    pattern = re.compile(
        r"^\s*ODE_StepSize\s*=\s*"
        r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
    )
    with world_file.open("r", encoding="utf-8", errors="replace") as stream:
        for line in stream:
            match = pattern.match(line.split("#", 1)[0])
            if match:
                value = float(match.group(1))
                if value <= 0.0:
                    fail(f"ODE_StepSize must be positive in {world_file}")
                return value
    fail(f"Cannot find ODE_StepSize in {world_file}")


def discover_steps(folder: Path, pattern: re.Pattern[str]) -> dict[int, Path]:
    result: dict[int, Path] = {}
    if not folder.is_dir():
        return result
    for path in folder.iterdir():
        if not path.is_file():
            continue
        match = pattern.match(path.name)
        if match:
            result[int(match.group(1))] = path
    return dict(sorted(result.items()))


def select_steps(
    files: dict[int, Path],
    start: int | None,
    end: int | None,
    stride: int,
    max_frames: int | None,
) -> list[tuple[int, Path]]:
    selected = [
        (step, path)
        for step, path in files.items()
        if (start is None or step >= start) and (end is None or step <= end)
    ]
    selected = selected[::stride]
    if max_frames is not None:
        selected = selected[:max_frames]
    return selected


def read_particle_state(path: Path) -> dict[str, array]:
    result = {
        "points": array("f"),
        "number": array("i"),
        "status": array("i"),
        "mass": array("f"),
        "radius": array("f"),
        "speed": array("f"),
        "angular_speed": array("f"),
        "velocity": array("f"),
        "angular_velocity": array("f"),
        "quaternion": array("f"),
    }
    with path.open("r", encoding="ascii") as stream:
        for line_number, line in enumerate(stream, 1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) < 17:
                fail(f"{path}:{line_number}: expected 17 columns, got {len(fields)}")
            number = int(fields[0])
            status = int(fields[1])
            mass = float(fields[2])
            radius = float(fields[3])
            position = tuple(float(value) for value in fields[4:7])
            velocity = tuple(float(value) for value in fields[7:10])
            quaternion = tuple(float(value) for value in fields[10:14])
            angular_velocity = tuple(float(value) for value in fields[14:17])

            result["number"].append(number)
            result["status"].append(status)
            result["mass"].append(mass)
            result["radius"].append(radius)
            result["points"].extend(position)
            result["velocity"].extend(velocity)
            result["quaternion"].extend(quaternion)
            result["angular_velocity"].extend(angular_velocity)
            result["speed"].append(math.sqrt(sum(value * value for value in velocity)))
            result["angular_speed"].append(
                math.sqrt(sum(value * value for value in angular_velocity))
            )
    if not result["number"]:
        fail(f"Particle state file is empty: {path}")
    return result


def write_particle_vtp(source: Path, destination: Path) -> int:
    data = read_particle_state(source)
    count = len(data["number"])
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="ascii", newline="\n") as stream:
        stream.write('<?xml version="1.0"?>\n')
        stream.write(
            '<VTKFile type="PolyData" version="1.0" byte_order="LittleEndian" '
            'header_type="UInt32">\n'
        )
        stream.write("  <PolyData>\n")
        stream.write(
            f'    <Piece NumberOfPoints="{count}" NumberOfVerts="{count}" '
            'NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">\n'
        )
        stream.write('      <PointData Scalars="radius" Vectors="velocity">\n')
        write_data_array(stream, "Int32", "number", 1, "i", data["number"], "        ")
        write_data_array(stream, "Int32", "status", 1, "i", data["status"], "        ")
        for name in ("mass", "radius", "speed", "angular_speed"):
            write_data_array(stream, "Float32", name, 1, "f", data[name], "        ")
        write_data_array(
            stream, "Float32", "velocity", 3, "f", data["velocity"], "        "
        )
        write_data_array(
            stream,
            "Float32",
            "angular_velocity_body",
            3,
            "f",
            data["angular_velocity"],
            "        ",
        )
        write_data_array(
            stream, "Float32", "quaternion_wxyz", 4, "f", data["quaternion"], "        "
        )
        stream.write("      </PointData>\n")
        stream.write("      <CellData/>\n")
        stream.write("      <Points>\n")
        write_data_array(stream, "Float32", None, 3, "f", data["points"], "        ")
        stream.write("      </Points>\n")
        stream.write("      <Verts>\n")
        write_data_array(stream, "Int32", "connectivity", 1, "i", range(count), "        ")
        write_data_array(stream, "Int32", "offsets", 1, "i", range(1, count + 1), "        ")
        stream.write("      </Verts>\n")
        stream.write("    </Piece>\n  </PolyData>\n</VTKFile>\n")
    return count


def read_bt_mesh(path: Path) -> BodyMesh:
    """Read the body surface mesh format used directly by the solver."""

    tokens = path.read_text(encoding="ascii").split()
    if len(tokens) < 2:
        fail(f"Invalid body mesh header: {path}")
    point_count = int(tokens[0])
    triangle_count = int(tokens[1])
    if point_count <= 0 or triangle_count <= 0:
        fail(f"Body mesh counts must be positive: {path}")

    expected_tokens = 2 + point_count * 3 + triangle_count * 3
    if len(tokens) < expected_tokens:
        fail(
            f"Incomplete body mesh {path}: expected at least "
            f"{expected_tokens} values, got {len(tokens)}"
        )

    cursor = 2
    coordinates = [float(value) for value in tokens[cursor : cursor + point_count * 3]]
    cursor += point_count * 3
    points = tuple(
        (coordinates[index], coordinates[index + 1], coordinates[index + 2])
        for index in range(0, len(coordinates), 3)
    )

    connectivity: list[int] = []
    for triangle in range(triangle_count):
        indices = [int(value) for value in tokens[cursor : cursor + 3]]
        cursor += 3
        if any(index < 1 or index > point_count for index in indices):
            fail(f"Body mesh triangle {triangle + 1} index out of range: {path}")
        connectivity.extend(index - 1 for index in indices)

    offsets = tuple(range(3, triangle_count * 3 + 1, 3))
    cell_types = (5,) * triangle_count  # VTK_TRIANGLE
    return BodyMesh(path.name, points, tuple(connectivity), offsets, cell_types)


def read_body_state(path: Path, expected_count: int) -> list[BodyState]:
    states: list[BodyState] = []
    with path.open("r", encoding="ascii") as stream:
        for line_number, line in enumerate(stream, 1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) < 14:
                fail(f"{path}:{line_number}: expected at least 14 columns, got {len(fields)}")
            states.append(
                BodyState(
                    tuple(float(value) for value in fields[1:4]),
                    tuple(float(value) for value in fields[10:14]),
                )
            )
    if len(states) != expected_count:
        fail(f"{path}: expected {expected_count} bodies, got {len(states)}")
    return states


def normalized_quaternion(
    quaternion: Sequence[float], conjugate: bool
) -> tuple[float, float, float, float]:
    norm = math.sqrt(sum(value * value for value in quaternion))
    if norm <= 1.0e-15:
        fail("Encountered a zero body quaternion")
    w, x, y, z = (value / norm for value in quaternion)
    return (w, -x, -y, -z) if conjugate else (w, x, y, z)


def rotation_matrix(quaternion: Sequence[float], conjugate: bool = False):
    w, x, y, z = normalized_quaternion(quaternion, conjugate)
    return (
        (1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - z * w), 2.0 * (x * z + y * w)),
        (2.0 * (x * y + z * w), 1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - x * w)),
        (2.0 * (x * z - y * w), 2.0 * (y * z + x * w), 1.0 - 2.0 * (x * x + y * y)),
    )


def rotate_vector(matrix, vector: Sequence[float]) -> tuple[float, float, float]:
    return tuple(sum(matrix[row][column] * vector[column] for column in range(3)) for row in range(3))


def transformed_point(
    matrix, point: Sequence[float], position: Sequence[float]
) -> tuple[float, float, float]:
    rotated = rotate_vector(matrix, point)
    return tuple(rotated[index] + position[index] for index in range(3))


def write_body_vtu(
    source: Path,
    destination: Path,
    meshes: Sequence[BodyMesh],
    conjugate_quaternions: bool,
) -> tuple[int, int]:
    states = read_body_state(source, len(meshes))
    points = array("f")
    connectivity = array("i")
    offsets = array("i")
    cell_types = array("B")
    point_body_ids = array("i")
    cell_body_ids = array("i")

    point_offset = 0
    connectivity_offset = 0
    for body_id, (mesh, state) in enumerate(zip(meshes, states)):
        matrix = rotation_matrix(state.quaternion, conjugate_quaternions)
        for point in mesh.points:
            points.extend(transformed_point(matrix, point, state.position))
        connectivity.extend(index + point_offset for index in mesh.connectivity)
        offsets.extend(value + connectivity_offset for value in mesh.offsets)
        cell_types.extend(mesh.cell_types)
        point_body_ids.extend([body_id] * len(mesh.points))
        cell_body_ids.extend([body_id] * len(mesh.cell_types))
        point_offset += len(mesh.points)
        connectivity_offset += len(mesh.connectivity)

    point_count = len(point_body_ids)
    cell_count = len(cell_body_ids)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="ascii", newline="\n") as stream:
        stream.write('<?xml version="1.0"?>\n')
        stream.write(
            '<VTKFile type="UnstructuredGrid" version="1.0" '
            'byte_order="LittleEndian" header_type="UInt32">\n'
        )
        stream.write("  <UnstructuredGrid>\n")
        stream.write(
            f'    <Piece NumberOfPoints="{point_count}" NumberOfCells="{cell_count}">\n'
        )
        stream.write('      <PointData Scalars="body_id">\n')
        write_data_array(stream, "Int32", "body_id", 1, "i", point_body_ids, "        ")
        stream.write("      </PointData>\n")
        stream.write('      <CellData Scalars="body_id">\n')
        write_data_array(stream, "Int32", "body_id", 1, "i", cell_body_ids, "        ")
        stream.write("      </CellData>\n")
        stream.write("      <Points>\n")
        write_data_array(stream, "Float32", None, 3, "f", points, "        ")
        stream.write("      </Points>\n")
        stream.write("      <Cells>\n")
        write_data_array(stream, "Int32", "connectivity", 1, "i", connectivity, "        ")
        write_data_array(stream, "Int32", "offsets", 1, "i", offsets, "        ")
        write_data_array(stream, "UInt8", "types", 1, "B", cell_types, "        ")
        stream.write("      </Cells>\n")
        stream.write("    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n")
    return point_count, cell_count


def write_pvd(
    destination: Path,
    entries: Sequence[tuple[int, Path]],
    time_step: float,
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write('<?xml version="1.0"?>\n')
        stream.write(
            '<VTKFile type="Collection" version="1.0" '
            'byte_order="LittleEndian">\n  <Collection>\n'
        )
        for step, path in entries:
            relative = path.relative_to(destination.parent).as_posix()
            stream.write(
                f'    <DataSet timestep="{step * time_step:.15g}" group="" '
                f'part="0" file="{relative}"/>\n'
            )
        stream.write("  </Collection>\n</VTKFile>\n")


def is_current(destination: Path, inputs: Sequence[Path], force: bool) -> bool:
    if force or not destination.is_file():
        return False
    destination_time = destination.stat().st_mtime
    return all(path.stat().st_mtime <= destination_time for path in inputs)


def load_body_meshes(body_folder: Path) -> tuple[list[BodyMesh], list[Path]]:
    """Load exactly the bodies declared by CouplePolyFileName.csv.

    The first CSV row is the authoritative body count. Each subsequent body
    row identifies the solver's .bt surface mesh. Legacy Data/DATA prefixes
    are intentionally ignored here; the matching file is loaded from the
    selected case's own BodySet directory.
    """

    csv_path = body_folder / "CouplePolyFileName.csv"
    if not csv_path.is_file():
        fail(f"Missing rigid-body list: {csv_path}")

    with csv_path.open("r", encoding="utf-8-sig", newline="") as stream:
        rows = csv.reader(stream)
        try:
            count_row = next(rows)
        except StopIteration:
            fail(f"Rigid-body list is empty: {csv_path}")
        try:
            body_count = int(count_row[0].strip())
        except (IndexError, ValueError):
            fail(f"Invalid rigid-body count in first row of {csv_path}")
        if body_count <= 0:
            fail(f"Rigid-body count must be positive in {csv_path}")

        try:
            header = next(rows)
        except StopIteration:
            fail(f"Missing header row in {csv_path}")
        header_names = [value.strip().lower() for value in header]
        if "filepath" not in header_names:
            fail(f"Cannot find filepath column in {csv_path}")
        filepath_column = header_names.index("filepath")

        mesh_paths: list[Path] = []
        for row_number, row in enumerate(rows, 3):
            if len(mesh_paths) == body_count:
                break
            if not row or not any(value.strip() for value in row):
                continue
            if filepath_column >= len(row) or not row[filepath_column].strip():
                fail(f"Missing filepath at {csv_path}:{row_number}")
            configured_path = Path(row[filepath_column].strip().replace("\\", "/"))
            mesh_path = body_folder / configured_path.name
            mesh_paths.append(mesh_path)

    if len(mesh_paths) != body_count:
        fail(
            f"{csv_path}: first row declares {body_count} bodies, "
            f"but only {len(mesh_paths)} body records were found"
        )
    missing = [str(path) for path in mesh_paths if not path.is_file()]
    if missing:
        fail("Missing body .bt mesh file(s):\n  " + "\n  ".join(missing))

    print(f"Bodies declared by CSV: {body_count}")
    print("Body meshes:")
    for index, path in enumerate(mesh_paths):
        print(f"  {index:2d}: {path.name}")
    meshes = [read_bt_mesh(path) for path in mesh_paths]
    return meshes, mesh_paths


def build_argument_parser() -> argparse.ArgumentParser:
    project_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description="Convert DEM/MBD .bt snapshots into particles.pvd and bodies.pvd."
    )
    parser.add_argument(
        "case_name",
        nargs="?",
        help="Case directory name under Data (for example DATA3)",
    )
    parser.add_argument(
        "--case",
        type=Path,
        help="Explicit case directory path (advanced override)",
    )
    parser.add_argument("--output", type=Path, help="Output directory (default: CASE/OutputFile/paraview)")
    parser.add_argument("--start-step", type=int)
    parser.add_argument("--end-step", type=int)
    parser.add_argument("--stride", type=int, default=1, help="Use every Nth available frame")
    parser.add_argument("--max-frames", type=int, help="Convert at most this many selected frames")
    parser.add_argument("--time-step", type=float, help="Override ODE_StepSize")
    parser.add_argument(
        "--conjugate-body-quaternions",
        action="store_true",
        help="Use conjugated body quaternions for legacy output compatibility",
    )
    parser.add_argument("--no-particles", action="store_true")
    parser.add_argument("--no-bodies", action="store_true")
    parser.add_argument(
        "--index-only",
        action="store_true",
        help="Rebuild PVD time indexes from existing VTK files without rewriting them",
    )
    parser.add_argument("--force", action="store_true", help="Rewrite files even if they are current")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_argument_parser().parse_args(argv)
    if args.stride <= 0:
        fail("--stride must be positive")
    if args.max_frames is not None and args.max_frames <= 0:
        fail("--max-frames must be positive")
    if args.no_particles and args.no_bodies:
        fail("Nothing to do: both --no-particles and --no-bodies were selected")

    project_root = Path(__file__).resolve().parents[1]
    config = read_case_config(Path(__file__).resolve().parent)
    configured_name = config.get("CASE_NAME", "DATA")
    if args.case is not None and args.case_name is not None:
        fail("Use either the positional CASE_NAME or --case, not both")
    case = (
        args.case
        if args.case is not None
        else case_directory_from_name(project_root, args.case_name or configured_name)
    ).resolve()
    output = (args.output or (case / "OutputFile" / "paraview")).resolve()
    world_file = case / "InputFile" / "SettingData" / "world.par"
    time_step = args.time_step if args.time_step is not None else read_step_size(world_file)
    if time_step <= 0.0:
        fail("Time step must be positive")
    converter_source = Path(__file__).resolve()
    print(f"Case:      {case}")
    print(f"Output:    {output}")
    print(f"Time step: {time_step:.15g} s")

    converted_any = False
    if not args.no_particles:
        particle_files = discover_steps(case / "OutputFile" / "state_particles", PARTICLE_PATTERN)
        selected = select_steps(
            particle_files, args.start_step, args.end_step, args.stride, args.max_frames
        )
        if not selected:
            print("Particle snapshots: none selected; particles.pvd was not generated.")
        else:
            print(f"Particle snapshots: {len(selected)}")
            pvd_entries: list[tuple[int, Path]] = []
            for index, (step, source) in enumerate(selected, 1):
                destination = output / "particles" / f"particles_{step:08d}.vtp"
                if args.index_only:
                    if not destination.is_file():
                        fail(f"Existing particle VTK file required by --index-only: {destination}")
                    action = "index"
                elif is_current(destination, (source, converter_source), args.force):
                    action = "reuse"
                else:
                    count = write_particle_vtp(source, destination)
                    action = f"write ({count} particles)"
                pvd_entries.append((step, destination))
                print(f"  particles {index:4d}/{len(selected)} step={step:8d}: {action}")
            write_pvd(output / "particles.pvd", pvd_entries, time_step)
            converted_any = True

    if not args.no_bodies:
        body_files = discover_steps(case / "OutputFile" / "state_bodys", BODY_PATTERN)
        selected = select_steps(
            body_files, args.start_step, args.end_step, args.stride, args.max_frames
        )
        if not selected:
            print("Body snapshots: none selected; bodies.pvd was not generated.")
        else:
            body_folder = case / "InputFile" / "BodySet"
            meshes, template_paths = load_body_meshes(body_folder)
            read_body_state(selected[0][1], len(meshes))
            conjugate = args.conjugate_body_quaternions
            print(
                "Body quaternion mode: "
                + ("conjugated wxyz" if conjugate else "direct wxyz")
            )
            print(f"Body snapshots: {len(selected)}")
            pvd_entries = []
            for index, (step, source) in enumerate(selected, 1):
                destination = output / "bodies" / f"bodies_{step:08d}.vtu"
                inputs = (source, converter_source, *template_paths)
                if args.index_only:
                    if not destination.is_file():
                        fail(f"Existing body VTK file required by --index-only: {destination}")
                    action = "index"
                elif is_current(destination, inputs, args.force):
                    action = "reuse"
                else:
                    point_count, cell_count = write_body_vtu(
                        source, destination, meshes, conjugate
                    )
                    action = f"write ({point_count} points, {cell_count} cells)"
                pvd_entries.append((step, destination))
                print(f"  bodies    {index:4d}/{len(selected)} step={step:8d}: {action}")
            write_pvd(output / "bodies.pvd", pvd_entries, time_step)
            converted_any = True

    if not converted_any:
        fail("No matching snapshots were found. Enable the corresponding world.par outputs first.")
    print("Done. Open the generated .pvd file(s) in ParaView.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
