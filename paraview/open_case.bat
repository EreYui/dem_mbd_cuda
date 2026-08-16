@echo off
setlocal EnableDelayedExpansion

set "SCRIPT_DIR=%~dp0"
set "PROJECT_ROOT=%SCRIPT_DIR%.."
set "CASE_NAME=DATA"
set "CONFIG_PARAVIEW_EXE="

if exist "%SCRIPT_DIR%case.dat" (
    for /f "usebackq tokens=1,* delims==" %%A in ("%SCRIPT_DIR%case.dat") do (
        if /i "%%A"=="CASE_NAME" set "CASE_NAME=%%B"
        if /i "%%A"=="PARAVIEW_EXE" set "CONFIG_PARAVIEW_EXE=%%B"
    )
)
if not "%~1"=="" (
    set "FIRST_ARGUMENT=%~1"
    if not "!FIRST_ARGUMENT:~0,1!"=="-" set "CASE_NAME=%~1"
)
set "CASE_DIR=%PROJECT_ROOT%\Data\%CASE_NAME%"

echo Selected case: %CASE_NAME%
echo Case directory: %CASE_DIR%

py -3 -c "import sys" >nul 2>nul
if not errorlevel 1 (
    py -3 "%SCRIPT_DIR%convert_case.py" %*
) else (
    python -c "import sys" >nul 2>nul
    if errorlevel 1 (
        echo ERROR: Python 3 was not found. Install Python 3 or add it to PATH.
        pause
        exit /b 1
    )
    python "%SCRIPT_DIR%convert_case.py" %*
)
if errorlevel 1 (
    echo.
    echo Conversion failed.
    pause
    exit /b 1
)

set "PV="
if defined PARAVIEW_EXE if exist "%PARAVIEW_EXE%" set "PV=%PARAVIEW_EXE%"
if not defined PV if defined CONFIG_PARAVIEW_EXE if exist "%CONFIG_PARAVIEW_EXE%" set "PV=%CONFIG_PARAVIEW_EXE%"
if not defined PV (
    for /f "delims=" %%P in ('where paraview.exe 2^>nul') do if not defined PV set "PV=%%P"
)
if not defined PV (
    for /d %%D in ("%ProgramFiles%\ParaView*") do (
        if exist "%%~fD\bin\paraview.exe" set "PV=%%~fD\bin\paraview.exe"
    )
)

if not defined PV (
    echo.
    echo Conversion completed, but paraview.exe was not found.
    echo Open these two files manually:
    echo   %CASE_DIR%\OutputFile\paraview\particles.pvd
    echo   %CASE_DIR%\OutputFile\paraview\bodies.pvd
    echo.
    echo Or set PARAVIEW_EXE to the full path of paraview.exe and run this file again.
    pause
    exit /b 2
)

echo Starting ParaView: %PV%
set "DEM_MBD_PARAVIEW_PROJECT_ROOT=%PROJECT_ROOT%"
set "DEM_MBD_PARAVIEW_CASE_DIR=%CASE_DIR%"
start "ParaView DEM-MBD" "%PV%" --script="%SCRIPT_DIR%open_case.py"
