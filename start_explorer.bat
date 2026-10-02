@echo off
rem Double-click to start the CML data explorer in your browser (stop it by closing this window).
rem Uses the miniconda Python if it exists (it has netCDF4, xarray, zarr, streamlit, plotly), else "python".
cd /d "%~dp0"
set "PY=%USERPROFILE%\miniconda3\python.exe"
if not exist "%PY%" set "PY=python"
"%PY%" -m cmlpipe explore %*
if errorlevel 1 pause
