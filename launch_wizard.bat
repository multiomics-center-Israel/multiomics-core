@echo off
:: ============================================================
:: multiomics-core - Launch Config Wizard
:: Double-click this file to open the wizard in your browser.
:: Bypasses renv activation for fast startup (~3s vs ~15s).
:: ============================================================
cd /d "%~dp0"
title Multiomics Pipeline - Loading...

:: If the bundled portable pandoc is present, point rmarkdown at it.
:: This avoids "error 127" when the install path contains spaces, since
:: RSTUDIO_PANDOC gets propagated verbatim and rmarkdown quotes it properly.
if exist "%~dp0tools\pandoc\pandoc.exe" (
    set "RSTUDIO_PANDOC=%~dp0tools\pandoc"
)

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
