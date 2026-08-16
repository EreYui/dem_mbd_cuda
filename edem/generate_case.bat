@echo off
setlocal
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0generate_case.ps1" %*
exit /b %errorlevel%
