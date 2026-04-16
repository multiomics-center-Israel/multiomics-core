@echo off
:: ============================================================
:: multiomics-core - Launch Config Wizard
:: Double-click this file to open the wizard in your browser.
:: Bypasses renv activation for fast startup (~3s vs ~15s).
:: ============================================================
cd /d "%~dp0"
title Multiomics Pipeline - Loading...
echo.
echo   Starting wizard - your browser will open shortly...
echo.
Rscript --vanilla run.R --wizard
if errorlevel 1 (
    echo.
    echo Something went wrong. Make sure R is installed and renv packages are restored.
    echo Run install.bat first if this is a fresh setup.
    pause
)
