@echo off
:: ============================================================
:: multiomics-core - Docker Launcher (Windows)
::
:: What this does:
::   1. Checks Docker Desktop is running
::   2. Resolves the image in this order:
::        a. Image already in Docker -> use it
::        b. multiomics-core_<TAG>.tar in the same folder -> docker load
::           (fastest offline route — pre-built on another machine)
::        c. Pull from ghcr.io -> needs the package to be Public
::        d. Build locally -> ~30 min, last resort
::   3. Starts the wizard on http://localhost:8080
::   4. Mounts your data + outputs folders so reports persist on Windows
::
:: Prerequisites:
::   - Docker Desktop for Windows installed and RUNNING
::     https://www.docker.com/products/docker-desktop
::   - WSL2 enabled (Docker Desktop installer handles this)
::
:: First run timing:
::   - .tar already next to this script:   ~5 min load + 10 sec start
::   - ghcr.io image publicly pullable:    ~5 min pull + 10 sec start
::   - Neither, but internet available:    ~30 min build + 10 sec start
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

:: Step 2: Resolve the image. Order of preference:
::   (a) image already in Docker     -> use it
::   (b) `.tar` file next to script  -> docker load it (offline transfer flow)
::   (c) ghcr.io pull                -> needs package to be public
::   (d) local docker build          -> last resort, ~30 min
docker image inspect "%IMAGE_LOCAL%" >nul 2>&1
if not errorlevel 1 (
    echo [OK] Using existing local image %IMAGE_LOCAL%.
    set "RUN_IMAGE=%IMAGE_LOCAL%"
    goto :run
)

:: (b) Look for a bundled .tar next to this script
set "TAR_FILE=%~dp0multiomics-core_%MULTIOMICS_IMAGE_TAG%.tar"
if exist "%TAR_FILE%" (
    echo Loading image from %TAR_FILE% ...
    docker load -i "%TAR_FILE%"
    if errorlevel 1 (
        echo [WARN] docker load failed. Will try ghcr.io / build next.
    ) else (
        echo [OK] Loaded image from .tar.
        docker image inspect "%IMAGE_LOCAL%" >nul 2>&1
        if not errorlevel 1 ( set "RUN_IMAGE=%IMAGE_LOCAL%" & goto :run )
        docker image inspect "%IMAGE_REMOTE%" >nul 2>&1
        if not errorlevel 1 ( set "RUN_IMAGE=%IMAGE_REMOTE%" & goto :run )
    )
)

:: (c) Pull from ghcr.io
echo Pulling prebuilt image from ghcr.io ...
echo   %IMAGE_REMOTE%
docker pull "%IMAGE_REMOTE%" 2>&1
if errorlevel 1 (
    echo.
    echo [WARN] Pull from ghcr.io failed. Falling back to local build.
    echo        Common reasons: no internet, package set to private, no .tar
    echo        bundled next to this script.
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

:run
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
