@echo off
:: ============================================================
:: multiomics-core — Launch Config Wizard
:: Double-click this file to open the wizard in your browser.
:: ============================================================
cd /d "%~dp0"
Rscript run.R --wizard
if errorlevel 1 (
    echo.
    echo Something went wrong. Make sure R is installed and renv packages are restored.
    echo Run install.bat first if this is a fresh setup.
    pause
)
