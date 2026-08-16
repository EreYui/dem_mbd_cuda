@echo off
setlocal

set "SCRIPT_DIR=%~dp0"
set "PROJECT_ROOT=%SCRIPT_DIR%.."

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

echo.
echo Conversion completed.
echo See the Output path printed by convert_case.py above.
pause
