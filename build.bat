@echo off
setlocal EnableExtensions EnableDelayedExpansion

for %%I in ("%~dp0.") do set "PROJECT_ROOT=%%~fI"
set "BUILD_DIR=%PROJECT_ROOT%\build-cuda"

rem Normalize duplicate Path/PATH entries injected by some launchers; MSBuild
rem rejects an environment containing both spellings.
set "SAVED_BUILD_PATH=%PATH%"
set "Path="
set "PATH=%SAVED_BUILD_PATH%"

where cmake >nul 2>nul || (
    echo [ERROR] CMake was not found on PATH.
    exit /b 1
)
where nvcc >nul 2>nul || (
    echo [ERROR] nvcc was not found on PATH. Install a supported CUDA Toolkit.
    exit /b 1
)
where nvidia-smi >nul 2>nul || (
    echo [ERROR] nvidia-smi was not found. Install an NVIDIA display driver.
    exit /b 1
)

where cl >nul 2>nul
if errorlevel 1 (
    set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
    if not exist "!VSWHERE!" (
        echo [ERROR] Visual Studio Installer vswhere.exe was not found.
        exit /b 1
    )
    for /f "usebackq tokens=*" %%I in (`"!VSWHERE!" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "VSROOT=%%I"
    if not defined VSROOT (
        echo [ERROR] Visual Studio 2022 C++ x64 tools were not found.
        exit /b 1
    )
    call "!VSROOT!\Common7\Tools\VsDevCmd.bat" -arch=x64 -host_arch=x64 || exit /b 1
)

echo [INFO] Toolchain
cmake --version | findstr /b /c:"cmake version"
nvcc --version | findstr /c:"release"
nvidia-smi --query-gpu=name,driver_version,compute_cap,memory.total --format=csv,noheader

cmake -S "%PROJECT_ROOT%" -B "%BUILD_DIR%" -G "Visual Studio 17 2022" -A x64 -DBUILD_TESTING=ON || exit /b 1
cmake --build "%BUILD_DIR%" --config Release --parallel || exit /b 1
ctest --test-dir "%BUILD_DIR%" -C Release --output-on-failure || exit /b 1
"%BUILD_DIR%\Release\dem_mbd.exe" --cuda-info || exit /b 1

echo [SUCCESS] Release executable: %BUILD_DIR%\Release\dem_mbd.exe
exit /b 0
