# ERA5 to LUCIE T30 Regridding Pipeline

Regrid ERA5 reanalysis data onto the SPEEDY T30 Gaussian grid (96x48) and
convert to per-hour HDF5 files matching the LUCIE training data format.

This pipeline reproduces the exact preprocessing used to create LUCIE's
training data, as described in:

> Arcomano, T., et al. (2022). "A Machine Learning-Based Global Atmospheric
> Forecast Model." *Geophysical Research Letters*, 49(9), e2021GL097268.

## Prerequisites

- **NCAR Derecho** HPC access with a compute allocation (e.g., `YOUR_ALLOCATION`)
- **Conda environment** with: `numpy`, `h5py`, `netCDF4`, `xarray`
- **Fortran compiler** via `ncarcompilers` module (uses `ftn` Cray wrapper)
- **CMake** for building the bspline-fortran dependency

## Quick Start

```bash
# 1. Clone with submodules
git clone --recurse-submodules <REPO_URL>
cd era5-lucie-regrid

# 2. Build the bspline library
module load ncarenv intel ncarcompilers cray-mpich netcdf cmake
cd bspline-fortran && mkdir -p build && cd build
cmake .. -DCMAKE_Fortran_COMPILER=ftn && make
cd ../..

# 3. Build the regridder
cd fortran && make && cd ..

# 4. Set up auxiliary files (symlinks to Troy's data — see below)
mkdir -p aux
ln -sf /glade/derecho/scratch/mdarman/ERA5_hr_haiwen/ERA5_T30_1950_TO_2021/speedy_average_surfacep/year_average_logp.nc aux/
ln -sf /glade/derecho/scratch/mdarman/ERA5_hr_haiwen/ERA5_T30_1950_TO_2021/orography_data/speedy_orography.nc aux/
ln -sf /glade/derecho/scratch/mdarman/ERA5_hr_haiwen/ERA5_T30_1950_TO_2021/orography_data/era_orography_jan01_1980.nc aux/

# 5. Generate data for one month
qsub -v YEAR=2022,MONTH=1,DO_H5=1 pbs/job_run_month.pbs
```

## Two Paths to Production

### Path A: Fast extract (recommended for 1950-2021)

Troy's bit-perfect yearly files already exist for 1950-2021. Extract per-hour
h5 directly — takes ~30 seconds per month:

```bash
qsub -v YEAR=2014 pbs/job_extract_year.pbs
```

### Path B: Full regridder (for years not in Troy's archive)

For new years (2022+), run the full pipeline: download ERA5 from RDA, regrid,
convert to h5:

```bash
qsub -v YEAR=2022,MONTH=1,DO_H5=1 pbs/job_run_month.pbs
```

See [RECIPE.md](RECIPE.md) for full details, verification results, and
troubleshooting.

## Configuration

All scripts support environment variable overrides so you can adapt to your
own paths without editing the scripts:

| Variable | Description | Default |
|----------|-------------|---------|
| `LUCIE_WORKDIR` | Root directory for inputs/outputs/aux | `/.../gen_lucie_3d` |
| `LUCIE_PYTHON` | Path to Python with numpy/h5py/netCDF4 | `/.../jax/bin/python` |
| `LUCIE_TROY_AUX` | Troy's yearly file archive | `/.../ERA5_T30_1950_TO_2021` |
| `LUCIE_TROY_YEARLY` | Troy's yearly NC files (Path A) | `/.../ERA5_T30_1950_TO_2021` |
| `LUCIE_REF_BASE` | LUCIE reference h5 files (for sst/tisr) | `/.../LUCIE_3D/1_step_1hr_h5df_test` |
| `LUCIE_SST_DIR` | SST boundary condition h5 files | `/.../sst_2014-2023` |
| `LUCIE_TISR_DIR` | TISR boundary condition h5 files | `/.../tisr_2015-2023` |

PBS files also need manual edits for:
- `#PBS -A <ALLOCATION>` — your compute allocation
- `#PBS -o <PATH>` — log output directory (PBS resolves this at parse time)

The Fortran source (`fortran/regrid_era.f90`) has a `base_dir` variable near
line 10 — change it to point to your working directory.

## Auxiliary Files

Three auxiliary NetCDF files are required in `aux/`:

| File | Description | Source |
|------|-------------|--------|
| `year_average_logp.nc` | Annual-mean log(surface pressure) on T30 grid | Troy's archive |
| `speedy_orography.nc` | SPEEDY model orography on T30 grid | Troy's archive |
| `era_orography_jan01_1980.nc` | ERA5 geopotential (orography) | Troy's archive |

**Critical**: `year_average_logp.nc` MUST be Troy's actual file. Reconstructing
it gives slightly different values that propagate to ~1 K temperature errors.

## Documentation

- [RECIPE.md](RECIPE.md) — Full production recipe with verification results
- [docs/horizontal_interpolation.md](docs/horizontal_interpolation.md) — How the 2D bspline regridding works
- [docs/vertical_interpolation.md](docs/vertical_interpolation.md) — How the sigma vertical interpolation works

## Repository Structure

```
era5-lucie-regrid/
├── README.md
├── RECIPE.md              # Production recipe
├── fortran/               # Fortran regridder source
├── bspline-fortran/       # B-spline library (git submodule)
├── scripts/
│   ├── download/          # CDS API download scripts (alternative to RDA)
│   ├── build_inputs_from_rda.py
│   ├── extract_h5_from_yearly.py
│   ├── monthly_to_lucie_h5.py
│   ├── run_month.sh
│   ├── run_range.sh
│   └── verify_against_reference.py
├── pbs/                   # PBS job submission scripts
└── docs/                  # Technical documentation
```
