#!/usr/bin/env python3
"""Generate an EDEM 2024 deck from one dem_mbd_omp case directory.

This script is intentionally executed with Altair's bundled Python 3.8 by
generate_case.ps1.  It uses the EDEMpy wheel shipped with EDEM 2024 and starts
from Altair's Linear_spring example so that the resulting deck has a valid
Linear Spring physics configuration.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import sys
from pathlib import Path

import numpy as np


EDEM_DEFAULT = Path(r"C:\Program Files\Altair\2024\EDEM")
TEMPLATE_RELATIVE = Path("examples") / "Linear_spring" / "Linear_spring"
DECK_BASENAME = "DATA_cube_drop"
CUBE_SIZE = 0.15
CUBE_MULTISPHERE_DIVISIONS = 4


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--case", default="DATA")
    parser.add_argument("--edem-root", type=Path, default=EDEM_DEFAULT)
    return parser.parse_args()


def require_file(path: Path) -> Path:
    if not path.is_file():
        raise FileNotFoundError("Required file not found: {}".format(path))
    return path


def read_world(path: Path) -> dict[str, float]:
    result: dict[str, float] = {}
    assignment = re.compile(
        r"^\s*([A-Za-z][A-Za-z0-9_]*)\s*=\s*"
        r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"
    )
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for line in stream:
            match = assignment.match(line)
            if match:
                result[match.group(1)] = float(match.group(2))
    required = [
        "ODE_StepSize",
        "ODE_EndStep",
        "ODE_OutputInterval",
        "MECH_PP_mu",
        "MECH_PP_mu_R",
        "MECH_PP_epsN",
        "MECH_PP_kN",
        "MECH_PT_mu",
        "MECH_PT_mu_R",
        "MECH_PT_epsN",
        "MECH_PW_mu",
        "MECH_PW_mu_R",
        "MECH_PW_epsN",
        "CONTROL_GRAVITY_X",
        "CONTROL_GRAVITY_Y",
        "CONTROL_GRAVITY_Z",
    ]
    missing = [name for name in required if name not in result]
    if missing:
        raise ValueError(
            "world.par is missing required numeric parameters: {}".format(
                ", ".join(missing)
            )
        )
    return result


def read_particles(path: Path) -> dict[str, np.ndarray]:
    data = np.loadtxt(str(path), dtype=np.float64)
    if data.ndim != 2 or data.shape[1] < 17:
        raise ValueError(
            "Particles.bt must contain at least 17 numeric columns; got shape {}".format(
                data.shape
            )
        )
    return {
        "ids": data[:, 0].astype(np.int64),
        "mass": data[:, 2],
        "radius": data[:, 3],
        "positions": data[:, 4:7],
        "velocities": data[:, 7:10],
        "orientations": data[:, 10:14],
        "angular_velocities": data[:, 14:17],
    }


def read_body(path: Path) -> dict[str, object]:
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as stream:
        rows = list(csv.reader(stream))
    if len(rows) < 3:
        raise ValueError("Rigid-body CSV has no body row: {}".format(path))
    row = rows[2]
    if len(row) < 24:
        raise ValueError("Rigid-body CSV row has only {} columns".format(len(row)))
    return {
        "source_mesh": row[0],
        "mass": float(row[1]),
        "position": [float(value) for value in row[2:5]],
        "velocity": [float(value) for value in row[5:8]],
        "angular_velocity": [float(value) for value in row[8:11]],
        # The source labels these lambda_x...lambda_w, but the C++ code and
        # Quat2DCM implementation consume the four values as (w, x, y, z).
        "orientation_wxyz": [float(value) for value in row[11:15]],
        "inertia": [
            [float(value) for value in row[15:18]],
            [float(value) for value in row[18:21]],
            [float(value) for value in row[21:24]],
        ],
    }


def quat_to_matrix(q: list[float]) -> list[list[float]]:
    w, x, y, z = q
    norm = math.sqrt(w * w + x * x + y * y + z * z)
    if norm == 0.0:
        raise ValueError("Rigid-body quaternion has zero length")
    w, x, y, z = (value / norm for value in (w, x, y, z))
    return [
        [
            w * w + x * x - y * y - z * z,
            2.0 * (x * y + w * z),
            2.0 * (x * z - w * y),
        ],
        [
            2.0 * (x * y - w * z),
            w * w - x * x + y * y - z * z,
            2.0 * (y * z + w * x),
        ],
        [
            2.0 * (x * z + w * y),
            2.0 * (y * z - w * x),
            w * w - x * x - y * y + z * z,
        ],
    ]

def container_mesh() -> tuple[list[list[float]], list[list[int]]]:
    # Five zero-thickness surfaces: bottom and four vertical walls, open top.
    quads = [
        # bottom, z = 0
        ([0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.5, 0.5, 0.0], [0.0, 0.5, 0.0]),
        # x = 0
        ([0.0, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.5, 0.6], [0.0, 0.0, 0.6]),
        # x = 0.5
        ([0.5, 0.0, 0.0], [0.5, 0.0, 0.6], [0.5, 0.5, 0.6], [0.5, 0.5, 0.0]),
        # y = 0
        ([0.0, 0.0, 0.0], [0.0, 0.0, 0.6], [0.5, 0.0, 0.6], [0.5, 0.0, 0.0]),
        # y = 0.5
        ([0.0, 0.5, 0.0], [0.5, 0.5, 0.0], [0.5, 0.5, 0.6], [0.0, 0.5, 0.6]),
    ]
    coords: list[list[float]] = []
    triangles: list[list[int]] = []
    for quad in quads:
        offset = len(coords)
        coords.extend([list(point) for point in quad])
        triangles.extend(
            [[offset, offset + 1, offset + 2], [offset, offset + 2, offset + 3]]
        )
    return coords, triangles


def write_ascii_stl(path: Path, solid: str, coords, triangles) -> None:
    with path.open("w", encoding="ascii", newline="\n") as stream:
        stream.write("solid {}\n".format(solid))
        for triangle in triangles:
            p0 = np.asarray(coords[triangle[0]], dtype=float)
            p1 = np.asarray(coords[triangle[1]], dtype=float)
            p2 = np.asarray(coords[triangle[2]], dtype=float)
            normal = np.cross(p1 - p0, p2 - p0)
            length = float(np.linalg.norm(normal))
            if length:
                normal /= length
            stream.write(
                "  facet normal {:.12e} {:.12e} {:.12e}\n".format(*normal)
            )
            stream.write("    outer loop\n")
            for point in (p0, p1, p2):
                stream.write(
                    "      vertex {:.12e} {:.12e} {:.12e}\n".format(*point)
                )
            stream.write("    endloop\n  endfacet\n")
        stream.write("endsolid {}\n".format(solid))


def replace_bed_particles_h5(
    deck_path: Path, particle_type_name: str, particles: dict[str, np.ndarray]
) -> None:
    """Replace template particle arrays in 0.h5.

    EDEMpy 1.5.0 on Windows crashes in its single-particle setter for this
    EDEM 2024 example deck. The HDF5 timestep schema is public through EDEMpy,
    so replacing these fixed-size datasets is both faster and more stable.
    """
    import h5py

    data_path = deck_path.parent / (deck_path.stem + "_data") / "0.h5"
    with h5py.File(str(data_path), "r+") as h5:
        timestep_root = h5["TimestepData"]
        timestep_key = next(iter(timestep_root.keys()))
        particle_types = timestep_root[timestep_key]["ParticleTypes"]
        particle_group = None
        for key in particle_types.keys():
            candidate = particle_types[key]
            name = candidate.attrs.get("name", b"")
            if isinstance(name, bytes):
                name = name.decode("utf-8")
            if name == particle_type_name:
                particle_group = candidate
                break
        if particle_group is None:
            raise KeyError(
                "Particle type {!r} not found in {}".format(
                    particle_type_name, data_path
                )
            )

        count = len(particles["positions"])
        zeros7 = np.zeros((count, 7), dtype=np.float64)
        replacements = {
            "ids": np.arange(1, count + 1, dtype=np.int32),
            "position": np.asarray(particles["positions"], dtype=np.float64),
            "velocity": np.asarray(particles["velocities"], dtype=np.float64),
            "angular velocity": np.asarray(
                particles["angular_velocities"], dtype=np.float64
            ),
            "orientation": np.asarray(
                particles["orientations"], dtype=np.float64
            ),
            "creation time": np.zeros(count, dtype=np.float64),
            "scale": np.ones(count, dtype=np.float64),
            "external force torque": zeros7,
            "force torque": zeros7.copy(),
            "us_force torque": zeros7.copy(),
        }
        for dataset_name, values in replacements.items():
            dataset = particle_group[dataset_name]
            if dataset.shape != values.shape:
                raise ValueError(
                    "{} has EDEM shape {}, expected {}".format(
                        dataset_name, dataset.shape, values.shape
                    )
                )
            dataset[...] = values
        particle_group.attrs.modify("size", count)

        # Keep the prototype metadata in TimestepData aligned with CreatorData.
        # The Altair template initially contains 25 mm spheres; EDEMpy updates
        # CreatorData but does not rewrite the existing timestep prototype.
        creator_types = h5["CreatorData"]["0"]["ParticleTypes"]
        creator_group = None
        for key in creator_types.keys():
            candidate = creator_types[key]
            name = candidate.attrs.get("name", b"")
            if isinstance(name, bytes):
                name = name.decode("utf-8")
            if name == particle_type_name:
                creator_group = candidate
                break
        if creator_group is None:
            raise KeyError("CreatorData particle prototype not found")

        radius = float(np.median(particles["radius"]))
        surface_area = 4.0 * math.pi * radius**2
        for group in (creator_group, particle_group):
            group.attrs["bounding box dimensions"] = np.asarray(
                [2.0 * radius] * 3, dtype=np.float64
            )
            group.attrs["raw surface area"] = surface_area
            group.attrs["sphericity"] = 1.0
        for attr_name, attr_value in creator_group.attrs.items():
            if attr_name != "size":
                particle_group.attrs[attr_name] = attr_value
        particle_group.attrs["size"] = count
        if "spheres" in creator_group:
            sphere_values = creator_group["spheres"][()]
            particle_group["spheres"][...] = sphere_values

        # Synchronize the cube prototype too. This prevents stale zero-sized
        # bounding metadata after creating a cuboid through EDEMpy.
        for creator_key in creator_types.keys():
            source = creator_types[creator_key]
            source_name = source.attrs.get("name", b"")
            if isinstance(source_name, bytes):
                source_name = source_name.decode("utf-8")
            if source_name != "FallingCube":
                continue
            target = None
            for timestep_key in particle_types.keys():
                candidate = particle_types[timestep_key]
                candidate_name = candidate.attrs.get("name", b"")
                if isinstance(candidate_name, bytes):
                    candidate_name = candidate_name.decode("utf-8")
                if candidate_name == source_name:
                    target = candidate
                    break
            if target is None:
                raise KeyError("FallingCube timestep data not found")
            surface_area = 6.0 * CUBE_SIZE**2
            volume = CUBE_SIZE**3
            for group in (source, target):
                group.attrs["bounding box dimensions"] = np.asarray(
                    [CUBE_SIZE, CUBE_SIZE, CUBE_SIZE], dtype=np.float64
                )
                group.attrs["raw surface area"] = surface_area
                group.attrs["raw volume"] = volume
                group.attrs["sphericity"] = (
                    math.pi ** (1.0 / 3.0)
                    * (6.0 * volume) ** (2.0 / 3.0)
                    / surface_area
                )
            target_size = int(target.attrs["size"])
            for attr_name, attr_value in source.attrs.items():
                if attr_name != "size":
                    target.attrs[attr_name] = attr_value
            target.attrs["size"] = target_size
            if "spheres" in source and "spheres" in target:
                target["spheres"][...] = source["spheres"][()]
        h5.flush()


def ensure_linear_spring_pairs_h5(deck_path: Path) -> None:
    import h5py

    def ensure_rows(h5, dataset_path: str, pairs: list[tuple[bytes, bytes, float]]) -> None:
        dataset = h5[dataset_path]
        data = list(dataset[()])
        existing = {(bytes(row["name1"]), bytes(row["name2"])) for row in data}
        changed = False
        for name1, name2, value in pairs:
            if (name1, name2) not in existing and (name2, name1) not in existing:
                data.append((name1, name2, value))
                changed = True
        if not changed:
            return
        parent_path, dataset_name = dataset_path.rsplit("/", 1)
        parent = h5[parent_path]
        dtype = dataset.dtype
        del parent[dataset_name]
        parent.create_dataset(dataset_name, data=np.asarray(data, dtype=dtype))

    data_path = deck_path.parent / (deck_path.stem + "_data") / "0.h5"
    with h5py.File(str(data_path), "r+") as h5:
        ensure_rows(
            h5,
            "CreatorData/0/PhysicsModels/SurfaceSurface/Linear Spring/data",
            [
                (b"CubeMaterial", b"plastic", 0.5),
                (b"CubeMaterial", b"CubeMaterial", 0.5),
            ],
        )
        ensure_rows(
            h5,
            "CreatorData/0/PhysicsModels/SurfaceGeometry/Linear Spring/data",
            [
                (b"equip plastic", b"CubeMaterial", 0.5),
            ],
        )
        h5.flush()


def cube_multisphere_definition(size: float, divisions: int) -> tuple[list[list[float]], float]:
    """Return sphere centers and radius for a simple cube-like rigid clump."""
    if divisions < 2:
        raise ValueError("A cube multisphere approximation needs at least two divisions")
    spacing = size / divisions
    radius = spacing / 2.0
    first = -0.5 * size + radius
    values = [first + i * spacing for i in range(divisions)]
    centers = [[x, y, z] for x in values for y in values for z in values]
    return centers, radius


def copy_template(edem_root: Path, output: Path) -> Path:
    source_base = edem_root / TEMPLATE_RELATIVE
    require_file(source_base.with_suffix(".dem"))
    output.mkdir(parents=True, exist_ok=True)
    target_base = output / DECK_BASENAME
    for generated_name in (DECK_BASENAME + ".efd", "edem_run.log"):
        generated_path = output / generated_name
        if generated_path.exists():
            generated_path.unlink()
    for suffix in (".dem", ".dfg", ".ess", ".ptf"):
        source = source_base.with_suffix(suffix)
        target = target_base.with_suffix(suffix)
        if source.exists():
            shutil.copy2(str(source), str(target))
    source_data = source_base.parent / (source_base.name + "_data")
    target_data = output / (DECK_BASENAME + "_data")
    if target_data.exists():
        shutil.rmtree(str(target_data))
    shutil.copytree(str(source_data), str(target_data))
    return target_base.with_suffix(".dem")


def create_deck(deck_path: Path, world, particles, body) -> None:
    import edempy

    def progress(message: str) -> None:
        print("[EDEM] {}".format(message), flush=True)

    particle_material = "plastic"
    equipment_material = "equip plastic"
    cube_material = "CubeMaterial"
    bed_type = "Ball"
    cube_type = "FallingCube"

    progress("opening template deck")
    with edempy.Deck(
        str(deck_path),
        mode="w",
        updateDeckVersion=True,
        warningBehavior="ignore",
    ) as deck:
        # Move all template geometry outside the active domain. EDEMpy 1.5 can
        # create but cannot delete geometry, so this keeps a valid official
        # template while making its original geometry inert.
        progress("moving template geometry")
        template_geometry_names = list(
            deck.creatorData[deck.numTimesteps - 1].geometryNames
        )
        template_motions = {
            name: list(
                deck.creatorData[deck.numTimesteps - 1]
                .geometry[name]
                .kinematicNames
            )
            for name in template_geometry_names
        }
        for geometry_name in template_geometry_names:
            deck.setGeometryTranslation(geometry_name, [100.0, 100.0, 100.0])
            for motion_name in template_motions[geometry_name]:
                deck.deleteKinematic(geometry_name, motion_name)

        mass = float(np.median(particles["mass"]))
        radius = float(np.median(particles["radius"]))
        density = mass / ((4.0 / 3.0) * math.pi * radius**3)
        particle_volume = (4.0 / 3.0) * math.pi * radius**3
        particle_moi = (2.0 / 5.0) * mass * radius**2

        # The shear modulus is selected so that EDEM's Linear Spring formula,
        # at the template characteristic velocity of 0.5 m/s, reproduces the
        # requested bed-particle normal stiffness as closely as possible.
        poisson = 0.25
        target_kn = world["MECH_PP_kN"]
        equivalent_radius = radius / 2.0
        equivalent_mass = mass / 2.0
        characteristic_velocity = 0.5
        equivalent_young = (
            target_kn**5
            / (
                (16.0 / 15.0) ** 4
                * equivalent_radius**2
                * equivalent_mass
                * characteristic_velocity**2
            )
        ) ** 0.25
        particle_shear = equivalent_young * (1.0 - poisson)

        progress("configuring materials and bed particle prototype")
        deck.setMaterialPoisson(particle_material, poisson)
        deck.setMaterialShearModulus(particle_material, particle_shear)
        deck.setMaterialDensity(particle_material, density)
        deck.setMaterialPoisson(equipment_material, poisson)
        deck.setMaterialShearModulus(equipment_material, particle_shear)
        deck.setMaterialDensity(equipment_material, 7800.0)

        deck.setMultisphereSpherePosition(bed_type, 0, [0.0, 0.0, 0.0])
        deck.setMultispherePhysicalRadius(bed_type, 0, radius)
        deck.setMultisphereContactRadius(bed_type, 0, radius)
        deck.setAutoCalculateParticleProperties(bed_type, False)
        deck.setParticleTypeMass(bed_type, mass)
        deck.setParticleTypeVolume(bed_type, particle_volume)
        deck.setParticleTypeMoI(
            bed_type, [particle_moi, particle_moi, particle_moi]
        )

        deck.setGravity(
            [
                world["CONTROL_GRAVITY_X"],
                world["CONTROL_GRAVITY_Y"],
                world["CONTROL_GRAVITY_Z"],
            ]
        )
        deck.setDomainMin([-0.05, -0.05, -0.02])
        deck.setDomainMax([0.55, 0.55, 0.70])

        progress("creating container geometry")
        coords, triangles = container_mesh()
        deck.createGeometry(
            "BedContainer", equipment_material, coords, triangles, "Physical"
        )

        # Grow the official EDEM datasets to the required size through EDEMpy.
        # Directly resizing/recreating these datasets is not accepted by the
        # native EDEM reader. The following in-place state replacement is fast
        # once EDEM has allocated the correct layout.
        progress("allocating bed particle records through EDEMpy")
        existing_count = int(
            deck.timestep[deck.numTimesteps - 1].particle[bed_type].numParticles
        )
        batch_size = 5000
        for start in range(existing_count, len(particles["positions"]), batch_size):
            end = min(start + batch_size, len(particles["positions"]))
            progress("allocating particles {}..{}".format(start, end - 1))
            deck.createParticles(
                bed_type,
                particles["positions"][start:end].tolist(),
                velocities=particles["velocities"][start:end].tolist(),
                ang_velocities=particles["angular_velocities"][start:end].tolist(),
                orientations=particles["orientations"][start:end].tolist(),
            )

        progress("creating 6-DOF polyhedral cube")
        cube_density = body["mass"] / (CUBE_SIZE**3)
        deck.createMaterial(
            cube_material,
            poisson,
            particle_shear,
            cube_density,
            0.0,
            True,
        )
        deck.createPolyhedralParticleType(cube_material, cube_type, "cuboid")
        deck.setPolyhedralCuboidDefinition(cube_type, CUBE_SIZE, CUBE_SIZE, CUBE_SIZE)
        deck.setAutoCalculateParticleProperties(cube_type, False)
        deck.setParticleTypeMass(cube_type, body["mass"])
        deck.setParticleTypeVolume(cube_type, CUBE_SIZE**3)
        inertia = body["inertia"]
        deck.setParticleTypeMoI(
            cube_type, [inertia[0][0], inertia[1][1], inertia[2][2]]
        )
        deck.createParticle(
            cube_type,
            body["position"],
            velocity=body["velocity"],
            ang_velocity=body["angular_velocity"],
            orientation=body["orientation_wxyz"],
        )

        # Existing template interaction pairs.
        progress("configuring interaction coefficients")
        deck.setInteractionRestitution(
            particle_material, particle_material, world["MECH_PP_epsN"]
        )
        deck.setInteractionStaticFriction(
            particle_material, particle_material, world["MECH_PP_mu"]
        )
        deck.setInteractionRollingFriction(
            particle_material, particle_material, world["MECH_PP_mu_R"]
        )
        deck.setInteractionRestitution(
            particle_material, equipment_material, world["MECH_PW_epsN"]
        )
        deck.setInteractionStaticFriction(
            particle_material, equipment_material, world["MECH_PW_mu"]
        )
        deck.setInteractionRollingFriction(
            particle_material, equipment_material, world["MECH_PW_mu_R"]
        )

        # New cube pairs. The bed-cube pair corresponds to the original PT
        # parameters, while cube-container uses PW.
        deck.createInteraction(
            particle_material,
            cube_material,
            world["MECH_PT_epsN"],
            world["MECH_PT_mu"],
            world["MECH_PT_mu_R"],
        )
        deck.createInteraction(
            cube_material,
            equipment_material,
            world["MECH_PW_epsN"],
            world["MECH_PW_mu"],
            world["MECH_PW_mu_R"],
        )
        deck.createInteraction(
            cube_material,
            cube_material,
            world["MECH_PT_epsN"],
            world["MECH_PT_mu"],
            world["MECH_PT_mu_R"],
        )
        progress("saving deck")
        deck.save()
    progress("deck generation complete")
    progress("replacing template particle arrays with source bed")
    replace_bed_particles_h5(deck_path, bed_type, particles)
    ensure_linear_spring_pairs_h5(deck_path)
    progress("particle import complete")


def main() -> int:
    args = parse_args()
    repo = args.repo.resolve()
    case_root = repo / "Data" / args.case
    input_root = case_root / "InputFile"
    output = case_root / "edem"
    world_path = require_file(input_root / "SettingData" / "world.par")
    particles_path = require_file(input_root / "Particles" / "Particles.bt")
    body_path = require_file(
        input_root / "BodySet" / "CouplePolyFileName.csv"
    )

    world = read_world(world_path)
    particles = read_particles(particles_path)
    body = read_body(body_path)
    if not np.allclose(particles["mass"], particles["mass"][0]):
        raise ValueError("EDEM converter currently requires one particle mass")
    if not np.allclose(particles["radius"], particles["radius"][0]):
        raise ValueError("EDEM converter currently requires one particle radius")

    deck_path = copy_template(args.edem_root, output)
    create_deck(deck_path, world, particles, body)

    coords, triangles = container_mesh()
    write_ascii_stl(output / "bed_container.stl", "BedContainer", coords, triangles)

    total_time = world["ODE_StepSize"] * world["ODE_EndStep"]
    write_interval = world["ODE_StepSize"] * world["ODE_OutputInterval"]
    manifest = {
        "source_case": args.case,
        "deck": deck_path.name,
        "edem_version": "2024",
        "model": {
            "bed_particle_type": "Ball",
            "bed_particle_count": int(len(particles["positions"])),
            "bed_particle_radius_m": float(particles["radius"][0]),
            "bed_particle_mass_kg": float(particles["mass"][0]),
            "bed_particle_density_kg_m3": float(
                particles["mass"][0]
                / ((4.0 / 3.0) * math.pi * particles["radius"][0] ** 3)
            ),
            "cube_particle_type": "FallingCube",
            "cube_size_m": [CUBE_SIZE, CUBE_SIZE, CUBE_SIZE],
            "cube_representation": "EDEM polyhedral cuboid",
            "cube_mass_kg": body["mass"],
            "cube_inertia_kg_m2": body["inertia"],
            "cube_position_m": body["position"],
            "cube_orientation_wxyz": body["orientation_wxyz"],
            "container_inner_size_m": [0.5, 0.5],
            "gravity_m_s2": [
                world["CONTROL_GRAVITY_X"],
                world["CONTROL_GRAVITY_Y"],
                world["CONTROL_GRAVITY_Z"],
            ],
        },
        "simulation": {
            "time_step_s": world["ODE_StepSize"],
            "end_step": int(world["ODE_EndStep"]),
            "total_time_s": total_time,
            "write_interval_s": write_interval,
            "recommended_engine": "CUDA",
        },
        "contact_mapping": {
            "bed_bed": "MECH_PP_*",
            "bed_cube": "MECH_PT_*",
            "bed_or_cube_container": "MECH_PW_*",
            "edem_base_model": "Linear Spring",
            "characteristic_velocity_m_s": 0.5,
        },
        "known_differences": [
            "The original rigid cube is represented by an EDEM polyhedral cuboid. Contact-model material pairs are completed in the generated HDF5 deck to ensure that the cube-bed and cube-container interactions use the same nominal coefficients as the source case.",
            "EDEM Linear Spring derives stiffness from material properties and characteristic velocity; MECH_PP_kN is matched at 0.5 m/s for identical bed particles.",
            "MECH_*_mu_T and MECH_*_beta have no direct parameter in EDEM's built-in Linear Spring plus Standard Rolling Friction combination.",
        ],
    }
    with (output / "case_manifest.json").open(
        "w", encoding="utf-8", newline="\n"
    ) as stream:
        json.dump(manifest, stream, ensure_ascii=False, indent=2)
        stream.write("\n")

    print("Generated EDEM deck: {}".format(deck_path))
    print("Bed particles: {}".format(len(particles["positions"])))
    print(
        "Run settings: step={} s, total={} s, write={} s".format(
            world["ODE_StepSize"], total_time, write_interval
        )
    )
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print("EDEM case generation failed: {}".format(exc), file=sys.stderr)
        raise
