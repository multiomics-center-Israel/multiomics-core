# Running multiomics-core in Docker

The Docker image bundles **R 4.5.2**, **Bioconductor 3.22**, all
`renv.lock` packages, the wizard, the three example datasets
(proteomics / metabolomics / lipidomics) and pandoc. Once the image is
built (or pulled), the only thing the user needs locally is **Docker
Desktop**.

## When to use Docker vs the native installer

| Use Docker | Use native `install.bat` |
|---|---|
| Different R version than 4.5.2 on the machine | Already on R 4.5.2 + RTools 4.5 |
| Bioconductor / RTools install kept failing | Native install worked previously |
| Want a clean isolated env per project | Need maximum performance / no virtualization overhead |
| Sharing the pipeline across multiple machines | Single-user, single-machine |

For Ifat-style "I just want it to work on my laptop", Docker is the
shorter path past the R-version / Bioconductor / RTools maze.

## One-time install (Windows)

1. Install **Docker Desktop**:
   <https://www.docker.com/products/docker-desktop>
   - Requires Windows 10/11 64-bit and WSL2 (the installer enables WSL2
     automatically if missing).
   - Reboot after install.
2. Start Docker Desktop. Wait for the tray icon to say **"Engine
   running"** (about 60 seconds on first boot).

Verify in Command Prompt:
```cmd
docker --version
docker info
```
Both should print without errors.

## Build (one-time, ~30 min)

From the repo root:

```cmd
docker-run.bat
```

If the image isn't local, this script builds it from this folder and
then starts the wizard. The build downloads R 4.5.2, all
`renv.lock` packages, and pandoc — that's where the 30 minutes goes.
**You only pay this once per renv.lock change.**

## Run (every subsequent time, ~10 sec)

Same script:

```cmd
docker-run.bat
```

It starts the wizard at <http://localhost:8080>. Open that URL in any
browser. Press **Ctrl+C** in the terminal window to stop.

By default the script mounts two folders next to itself:

| Host folder              | Container path     | Purpose                                |
|--------------------------|--------------------|----------------------------------------|
| `<repo>\user_data\`      | `/app/data/user`   | Your input files (CSV/TSV/XLSX)        |
| `<repo>\outputs\`        | `/app/outputs`     | Generated reports + intermediate files |

Override by setting `USER_DATA_DIR` and/or `OUTPUTS_DIR` before
launching:

```cmd
set USER_DATA_DIR=C:\Users\Ifat\my_metabolomics_run\data
set OUTPUTS_DIR=C:\Users\Ifat\my_metabolomics_run\outputs
docker-run.bat
```

## Trying the example datasets

In the wizard:

1. Click **"Load Example (Metabolomics)"** (or Lipidomics, or Proteomics).
2. Scroll down. **Run Pipeline**.
3. When the run finishes, the HTML report opens automatically and a
   copy lives in `<repo>\outputs\` on Windows.

No file path setup required — the example data is baked into the image.

## Using your own data

1. Copy your input files into `<repo>\user_data\` on Windows.
   - e.g. `user_data\my_run\metab_data.csv`,
     `user_data\my_run\metadata.csv`
2. In the wizard, reference them with the in-container path
   `data/user/my_run/metab_data.csv`.
3. Run.

## Distributing the image without rebuilding

To save Ifat the 30-minute first build:

**On a machine where the image is already built:**
```cmd
docker save multiomics-core:latest -o multiomics-core.tar
```
That writes a ~3-5 GB `.tar`. Transfer it to her laptop (USB, cloud
drive, whatever).

**On her laptop:**
```cmd
docker load -i multiomics-core.tar
```
That loads the image into her local Docker. Then `docker-run.bat`
skips the build step and just starts.

For a more permanent solution, the image can be pushed to
[GitHub Container Registry](https://ghcr.io) — out of scope for this
file; ask if you want the GitHub Actions workflow.

## Troubleshooting

- **"Docker daemon not responding"** — Docker Desktop isn't running.
  Start it and wait for "Engine running".
- **"Port 8080 already in use"** — another process has 8080. Stop it or
  edit `docker-run.bat` to use a different host port
  (`-p 8081:8080` then go to <http://localhost:8081>).
- **Wizard loads but Run Pipeline hangs** — check `<repo>\outputs\` for
  partial files and the terminal output for stack traces. The
  container's `outputs/` directory is bind-mounted, so logs persist on
  the host.
