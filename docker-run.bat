@echo off
:: ============================================================
:: multiomics-core - Docker Launcher (Windows)
::
:: What this does:
::   1. Checks Docker Desktop is running
::   2. Pulls the prebuilt image from ghcr.io (~5 min, no R/Bioc build)
::      Falls back to a local build only if the pull fails (e.g. offline,
::      or the registry image hasn't been published yet for this branch).
::   3. Starts the wizard on http://localhost:8080
::   4. Mounts your data + outputs folders so reports persist on Windows
::
:: Prerequisites:
::   - Docker Desktop for Windows installed and RUNNING
::     https://www.docker.com/products/docker-desktop
::   - WSL2 enabled (Docker Desktop installer handles this)
::
:: First run on a clean machine:
::   - With internet:  ~5 min image pull, then ~10 sec to start.
::   - Without internet (or registry image missing): falls back to local
::     build, ~30 min.
:: Subsequent runs: ~10 sec.
:: ============================================================

setlocal
cd /d "%~dp0"

:: Image identity. Override by setting MULTIOMICS_IMAGE_TAG before launching
:: (e.g. set MULTIOMICS_IMAGE_TAG=v1.0 for a pinned release).
if not defined MULTIOMICS_IMAGE_TAG set "MULTIOMICS_IMAGE_TAG=Lipidomics_dor"
set "IMAGE_REMOTE=ghcr.io/multiomics-center-israel/multiomics-core:%MULTIOMICS_IMAGE_TAG%"
set "IMAGE_LOCAL=multiomics-core:%MULTIOMICS_IMAGE_TAG%"

echo.
echo ===========================================================
echo   multiomics-core - Docker Launcher
echo   Image tag: %MULTIOMICS_IMAGE_TAG%
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

:: Step 2: Try to pull the prebuilt image from ghcr.io first (fast path).
::         Fall back to local build if pull fails.
echo Pulling prebuilt image from ghcr.io ...
echo   %IMAGE_REMOTE%
docker pull "%IMAGE_REMOTE%" 2>&1
if errorlevel 1 (
    echo.
    echo [WARN] Pull from ghcr.io failed. Falling back to local build.
    echo        Reasons can include: no internet, registry image not yet
    echo        published for this branch, or repo set to private.
    echo.
    docker image inspect "%IMAGE_LOCAL%" >nul 2>&1
    if errorlevel 1 (
        echo Building local image (~30 minutes first time, cached afterwards) ...
        docker build -t "%IMAGE_LOCAL%" .
        if errorlevel 1 (
            echo.
            echo [ERROR] Build failed. Scroll up for the failing step.
            pause
            exit /b 1
        )
        echo [OK] Local image built.
    ) else (
        echo [OK] Using existing local image %IMAGE_LOCAL%.
    )
    set "RUN_IMAGE=%IMAGE_LOCAL%"
) else (
    echo [OK] Pulled %IMAGE_REMOTE%.
    set "RUN_IMAGE=%IMAGE_REMOTE%"
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
echo Starting wizard at http://localhost:8080
echo Open that URL in your browser. Press Ctrl+C in this window to stop.
echo.

:: Try to open the browser ~8 seconds after startup
start "" /b cmd /c "timeout /t 8 /nobreak >nul && start http://localhost:8080"

docker run --rm -it ^
    -p 8080:8080 ^
    -v "%USER_DATA_DIR%:/app/data/user" ^
    -v "%OUTPUTS_DIR%:/app/outputs" ^
    --name multiomics-core-wizard ^
    "%RUN_IMAGE%"

endlocal
