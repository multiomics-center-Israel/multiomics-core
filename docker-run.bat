@echo off
:: ============================================================
:: multiomics-core - Docker Launcher (Windows)
::
:: What this does:
::   1. Checks Docker Desktop is running
::   2. Pulls / builds the multiomics-core image if not present
::   3. Starts the wizard on http://localhost:8080
::   4. Mounts your data + outputs folders so reports persist on Windows
::
:: Prerequisites:
::   - Docker Desktop for Windows installed and RUNNING
::     https://www.docker.com/products/docker-desktop
::   - WSL2 enabled (Docker Desktop installer handles this)
::
:: First run:
::   - If the image isn't local, this script builds it from this folder.
::   - Build takes ~30 minutes (Bioconductor packages). Once. After that:
::     ~10 seconds to start the wizard.
:: ============================================================

setlocal
cd /d "%~dp0"

echo.
echo ===========================================================
echo   multiomics-core - Docker Launcher
echo ===========================================================
echo.

:: Step 1: Docker available?
where docker >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Docker not found on PATH.
    echo.
    echo Install Docker Desktop:
    echo   https://www.docker.com/products/docker-desktop
    echo.
    echo After install, start Docker Desktop and wait for "Engine running".
    echo Then run this script again.
    pause
    exit /b 1
)

:: Probe the daemon (returns non-zero if Docker Desktop isn't running)
docker info >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Docker daemon not responding.
    echo.
    echo Start Docker Desktop and wait for "Engine running" in the tray,
    echo then re-run this script.
    pause
    exit /b 1
)

echo [OK] Docker daemon is up.
echo.

:: Step 2: Image present? Build if not.
docker image inspect multiomics-core:latest >nul 2>&1
if errorlevel 1 (
    echo Image multiomics-core:latest not found locally. Building...
    echo This takes ~30 minutes the first time. Only happens once.
    echo.
    docker build -t multiomics-core:latest .
    if errorlevel 1 (
        echo [ERROR] Build failed. Scroll up for the failing step.
        pause
        exit /b 1
    )
    echo.
    echo [OK] Image built.
) else (
    echo [OK] Image multiomics-core:latest is present.
)
echo.

:: Step 3: Mount user data + outputs folders. Defaults to ./user_data and
:: ./outputs next to this script. Override by setting USER_DATA_DIR and/or
:: OUTPUTS_DIR before launching.
if not defined USER_DATA_DIR set "USER_DATA_DIR=%~dp0user_data"
if not defined OUTPUTS_DIR  set "OUTPUTS_DIR=%~dp0outputs"
if not exist "%USER_DATA_DIR%" mkdir "%USER_DATA_DIR%"
if not exist "%OUTPUTS_DIR%"  mkdir "%OUTPUTS_DIR%"

echo Mounting:
echo   Windows %USER_DATA_DIR%   -^> container /app/data/user
echo   Windows %OUTPUTS_DIR%    -^> container /app/outputs
echo.

:: Step 4: Run.
:: --rm        : remove container when the wizard exits
:: -p 8080:8080: forward wizard port to host
:: -v          : bind-mount data + outputs so reports persist on Windows
:: --name      : friendly name so re-runs don't pile up containers
echo Starting wizard at http://localhost:8080
echo Open that URL in your browser. Press Ctrl+C in this window to stop.
echo.
docker run --rm -it ^
    -p 8080:8080 ^
    -v "%USER_DATA_DIR%:/app/data/user" ^
    -v "%OUTPUTS_DIR%:/app/outputs" ^
    --name multiomics-core-wizard ^
    multiomics-core:latest

endlocal
