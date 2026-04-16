@echo off
:: ============================================================
:: multiomics-core — One-Time Setup
::
:: What this does:
::   1. Checks that R is installed
::   2. Installs all required R packages (via renv)
::   3. Creates a desktop shortcut for the wizard
::
:: Prerequisites:
::   - R 4.5+ installed (https://cran.r-project.org/bin/windows/base/)
::   - RTools 4.5 installed (https://cran.r-project.org/bin/windows/Rtools/)
::
:: ============================================================

echo.
echo ===========================================================
echo   multiomics-core — Setup
echo ===========================================================
echo.

cd /d "%~dp0"

:: Step 1: Check R is available
where Rscript >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Rscript not found on PATH.
    echo.
    echo Please install R first:
    echo   https://cran.r-project.org/bin/windows/base/
    echo.
    echo During installation, make sure to check:
    echo   "Add R to PATH" or "Save version number in registry"
    echo.
    echo After installing R, also install RTools:
    echo   https://cran.r-project.org/bin/windows/Rtools/
    echo.
    pause
    exit /b 1
)

echo [OK] R found:
Rscript --version 2>&1
echo.

:: Step 2: Restore renv packages
echo [2/3] Installing R packages (this may take 10-20 minutes on first run)...
echo.
Rscript -e "if (!requireNamespace('renv', quietly=TRUE)) install.packages('renv', repos='https://cloud.r-project.org'); renv::restore(prompt=FALSE)"

if errorlevel 1 (
    echo.
    echo [WARNING] Package installation had issues. Some packages may need RTools.
    echo Make sure RTools is installed: https://cran.r-project.org/bin/windows/Rtools/
    echo Then run this script again.
    echo.
    pause
    exit /b 1
)

echo.
echo [OK] All packages installed.
echo.

:: Step 3: Create desktop shortcut
echo [3/3] Creating desktop shortcut...
cscript //nologo "%~dp0create_shortcut.vbs"
echo.

echo ===========================================================
echo   Setup complete!
echo.
echo   You can now launch the wizard by:
echo     - Double-clicking the "Multiomics Pipeline" icon on your desktop
echo     - Or double-clicking launch_wizard.bat in this folder
echo ===========================================================
echo.
pause
