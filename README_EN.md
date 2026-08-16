# DEM–MBD CUDA

**English** | [中文](README.md)

A general-purpose, double-precision coupled Discrete Element Method (DEM) and Multibody Dynamics (MBD) simulator. NVIDIA CUDA is the only production compute backend. If no usable CUDA device is available, the program reports a clear error and exits; it never silently falls back to CPU simulation.

## Features

- `CONTROL_MULTIBODY_FLAG = 0`: pure DEM simulation.
- `CONTROL_MULTIBODY_FLAG = 1`: general DEM–MBD coupling, with body, triangle, and particle counts determined dynamically from input files.
- Each body independently selects force-driven, fixed, or prescribed motion through the trailing `state` column in `CouplePolyFileName.csv`.
- Particle SoA data, contact histories, sorted spatial grids, forces, and torques remain resident on the GPU.
- CUDA performs neighbor search; particle–particle, particle–wall, and particle–triangle contact; body-load reduction; and particle leapfrog integration.
- The CPU handles configuration and file I/O, the small general MBD system, and limited per-step body-state/load exchange.
- Existing parameter names, input layouts, output filenames, field order, and step numbering are preserved.
- Both `/` and `\` input separators are accepted. Case selectors such as `DATA` and `DATA1` remain supported.
- Output uses a two-frame `std::thread` pipeline. Disabled output types allocate no large body-by-particle buffers and trigger no corresponding full-state downloads.

See the [CUDA migration and validation report](docs/CUDA_MIGRATION_VALIDATION.md) for architecture, numerical validation, performance measurements, and known limitations.

## Requirements

- An NVIDIA GPU, a compatible display driver, and the CUDA Toolkit.
- CMake 3.24 or newer.
- C++17 and CUDA C++17 toolchains.
- Windows: Visual Studio 2022 x64 C++ tools.
- Ubuntu: Ubuntu 24.04 and a GCC/G++ version supported by the installed CUDA Toolkit.

The default architecture is `CMAKE_CUDA_ARCHITECTURES=native`; users may override it during configuration. The project does not enable `--use_fast_math`, and simulation data uses `double` precision by default.

## Build and Run on Windows

From the repository root:

```bat
build.bat
```

The script checks CMake, `nvcc`, the NVIDIA driver/GPU, and the Visual Studio 2022 x64 environment. It then creates `build-cuda`, builds Release, runs CTest, and probes the CUDA device.

```bat
.\build-cuda\Release\dem_mbd.exe --cuda-info
.\build-cuda\Release\dem_mbd.exe DATA
.\build-cuda\Release\dem_mbd.exe DATA1
.\build-cuda\Release\dem_mbd.exe F:\cases\my_case
.\build-cuda\Release\dem_mbd.exe F:\cases\my_case\InputFile\SettingData\world.par
```

## Build and Run on Ubuntu 24.04

```bash
bash build_ubuntu.sh
./build-cuda-linux/dem_mbd DATA
./build-cuda-linux/dem_mbd DATA1
```

To build, test, and immediately launch a case:

```bash
bash build_ubuntu.sh --run DATA1
```

Example with an explicit GPU architecture:

```bash
cmake -S . -B build-cuda-linux -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CUDA_ARCHITECTURES=120 -DBUILD_TESTING=ON
cmake --build build-cuda-linux --parallel
ctest --test-dir build-cuda-linux --output-on-failure
```

## Cases and Output Safety

Passing only `DATA1` resolves to `Data/DATA1`. A case directory or a direct path to `world.par` may also be supplied. Legacy paths such as `Data\DATA\...` inside parameter files are redirected to the selected case directory.

The simulator writes into that case's `OutputFile` directory using the existing naming convention. Copy inputs to a separate validation directory before testing if existing results must be preserved. Automated tests use `validation/cases` and dedicated short-run directories rather than production `Data/DATA*` outputs.

Large production case data and generated validation output are intentionally not versioned. Add or copy a case under `Data/DATA`, `Data/DATA1`, and so on before running it.

## Per-Body Motion State

`CouplePolyFileName.csv` accepts an optional 25th column named `state`, following the original 24 columns:

```csv
filepath,mass,x,y,z,vx,vy,vz,omega_x,omega_y,omega_z,lambda_x,lambda_y,lambda_z,lambda_w,I1,I2,I3,I4,I5,I6,I7,I8,I9,state
Data\DATA\InputFile\BodySet\plate1.bt,1,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0
```

- `state=0`: dynamic body whose motion is determined by contact, gravity, and other applied loads.
- `state=1`: fixed body. Its center of mass and orientation remain unchanged, and velocities supplied in the CSV are zeroed. Contact loads are still computed and may be written to output.
- `state=2`: prescribed body moving with constant `vx,vy,vz` and `omega_x,omega_y,omega_z` from the CSV. Contact loads are computed but do not alter the prescribed motion.

Legacy 24-column files remain valid and default to `state=0`. A non-integer value or an integer outside `0..2` produces an input error containing the absolute filename and body row.

## Tests

```bat
ctest --test-dir build-cuda -C Release --output-on-failure
```

The suite covers contact geometry, cross-platform path handling, CUDA device detection, free fall, one-particle wall contact, two-particle contact, one-particle triangle contact, all three body motion states, a short pure-DEM case, and a short dynamic two-body DEM–MBD case.

Output comparison utilities:

```bat
.\build-cuda\Release\output_compare.exe reference_output_dir cuda_output_dir
.\build-cuda\Release\state_diagnostics.exe particle_state_file body_state_file
```

## Diagnostics

Input failures report the input's purpose, resolved absolute path, and specific reason. CUDA kernel launches are checked for identifiable errors. To diagnose illegal memory access:

```bat
compute-sanitizer --tool memcheck --error-exitcode 99 ^
  .\build-cuda\Release\dem_mbd.exe validation\cases\triangle_contact
```

## Documentation

- [CUDA migration and validation report](docs/CUDA_MIGRATION_VALIDATION.md)
- [Changelog](CHANGELOG.md)
- [EDEM comparison notes](edem/README.md)
