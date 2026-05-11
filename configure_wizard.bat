@echo off
:: ============================================================
:: multiomics-core - Configure Wizard Defaults
::
:: Run this once during installation to customize the wizard
:: for a specific lab/center. Sets which modes are available,
:: default values, theme, etc.
:: ============================================================

cd /d "%~dp0"

echo.
echo ===========================================================
echo   multiomics-core - Wizard Configuration
echo ===========================================================
echo.
echo This will customize the wizard defaults for this installation.
echo You can re-run this anytime to change settings.
echo.

Rscript --vanilla configure_wizard.R

if errorlevel 1 (
    echo.
    echo Configuration failed. Make sure R is installed.
    pause
    exit /b 1
)

echo.
pause
