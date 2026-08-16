#!/usr/bin/env bash
set -Eeuo pipefail

project_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
build_dir="${BUILD_DIR:-$project_root/build-cuda-linux}"
jobs="${JOBS:-$(nproc 2>/dev/null || printf '1')}"
run_case=""

usage() {
    printf '%s\n' \
        "Usage: bash build_ubuntu.sh [--run CASE] [--build-dir PATH] [--jobs N]" \
        "Builds the CUDA-only Release target, runs CTest, and probes the GPU."
}

while (($#)); do
    case "$1" in
        --run) [[ $# -ge 2 ]] || { printf 'Missing value for --run\n' >&2; exit 2; }; run_case="$2"; shift 2 ;;
        --build-dir) [[ $# -ge 2 ]] || { printf 'Missing value for --build-dir\n' >&2; exit 2; }; build_dir="$2"; shift 2 ;;
        --jobs) [[ $# -ge 2 ]] || { printf 'Missing value for --jobs\n' >&2; exit 2; }; jobs="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) printf 'Unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ -r /etc/os-release ]]; then
    # shellcheck disable=SC1091
    source /etc/os-release
    if [[ "${ID:-}" != ubuntu || "${VERSION_ID:-}" != 24.04 ]]; then
        printf 'Warning: designed for Ubuntu 24.04; detected %s %s.\n' \
            "${ID:-unknown}" "${VERSION_ID:-unknown}" >&2
    fi
fi

for command_name in cmake nvcc nvidia-smi c++; do
    command -v "$command_name" >/dev/null 2>&1 || {
        printf 'Error: required command not found: %s\n' "$command_name" >&2
        exit 1
    }
done

printf 'Toolchain:\n'
cmake --version | head -n 1
nvcc --version | tail -n 1
nvidia-smi --query-gpu=name,driver_version,compute_cap,memory.total --format=csv,noheader

cmake -S "$project_root" -B "$build_dir" \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=ON
cmake --build "$build_dir" --parallel "$jobs"
ctest --test-dir "$build_dir" --output-on-failure
"$build_dir/dem_mbd" --cuda-info

printf 'Release executable: %s/dem_mbd\n' "$build_dir"
if [[ -n "$run_case" ]]; then
    cd -- "$project_root"
    "$build_dir/dem_mbd" "$run_case"
fi
