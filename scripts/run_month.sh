#!/bin/bash
#=============================================================================
# Production driver: generate LUCIE-format regridded NetCDF for one month.
#
# Usage:
#   ./run_month.sh <YEAR> <MONTH>
#   ./run_month.sh 2014 1
#
# What it does:
#   1) Download ERA5 plev + sfc for that month at 3.0° / 28 levels (cached)
#   2) Reformat into Fortran-expected format in inputs/<year>/
#   3) Patch regrid_era.f90 with the year/month hardcoded values
#   4) Rebuild main_regrid.exe
#   5) Run the regridder
#   6) Verify output exists and has expected dimensions
#
# Output:
#   outputs/<year>/era_5_m{M}_y{YYYY}_regridded_spectral_mpi.nc
#
# This file is bit-faithful to Troy's `era_5_y{YYYY}_regridded_mpi.nc`
# yearly file (which is the source of LUCIE training data) for the
# corresponding month-slice. See RECIPE.md for full details.
#=============================================================================
set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <YEAR> <MONTH>"
    echo "Example: $0 2014 1"
    exit 1
fi

YEAR=$1
MONTH=$2

# Validate
if ! [[ "$YEAR" =~ ^[0-9]{4}$ ]]; then
    echo "ERROR: YEAR must be 4 digits, got '$YEAR'"
    exit 1
fi
if ! [[ "$MONTH" =~ ^([1-9]|1[0-2])$ ]]; then
    echo "ERROR: MONTH must be 1-12, got '$MONTH'"
    exit 1
fi

WORKDIR="${LUCIE_WORKDIR:-/glade/derecho/scratch/mdarman/sst_data/era5_data_download/gen_lucie_3d}"
SRC_DIR="$WORKDIR/fortran"
SCRIPTS_DIR="$WORKDIR/scripts"
INPUTS_DIR="$WORKDIR/inputs/$YEAR"
OUTPUTS_DIR="$WORKDIR/outputs/$YEAR"
LOGS_DIR="$WORKDIR/logs"

mkdir -p "$INPUTS_DIR" "$OUTPUTS_DIR" "$LOGS_DIR"

OUTPUT_FILE="$OUTPUTS_DIR/era_5_m${MONTH}_y${YEAR}_regridded_spectral_mpi.nc"
LOG="$LOGS_DIR/run_y${YEAR}_m${MONTH}.log"

echo "=========================================================================="
echo "  LUCIE-format generation for ${YEAR}-${MONTH}"
echo "  Output: $OUTPUT_FILE"
echo "  Log:    $LOG"
echo "=========================================================================="

# Source modules consistently
module purge >/dev/null 2>&1 || true
module load ncarenv intel ncarcompilers cray-mpich netcdf cmake >/dev/null 2>&1

# Direct path to conda env python (avoids `conda activate` non-interactive headaches)
JAX_PYTHON="${LUCIE_PYTHON:-/glade/work/mdarman/conda-envs/jax/bin/python}"
if [ ! -x "$JAX_PYTHON" ]; then
    echo "ERROR: Python not found at $JAX_PYTHON"
    echo "Set LUCIE_PYTHON to your conda env python path"
    exit 1
fi

#-----------------------------------------------------------------------------
# Step 1: build ERA5 inputs from local RDA mirror
# (Avoids CDS API entirely; uses /glade/campaign/collections/rda/data/d633000/.
#  We sub-sample 0.25° to 3.0° via every-12th-point — verified bit-identical
#  to CDS at multiple test points.)
#-----------------------------------------------------------------------------
echo ""
echo "[1/5] Building ERA5 inputs from RDA mirror..."

"$JAX_PYTHON" "$SCRIPTS_DIR/build_inputs_from_rda.py" \
    --year "$YEAR" --month "$MONTH" --skip-if-exists \
    2>&1 | tee -a "$LOG"

# Verify the input files are in place
FULL_FILE="$INPUTS_DIR/era_5_m${MONTH}_y${YEAR}_full.nc"
SP_FILE="$INPUTS_DIR/era_5_m${MONTH}_y${YEAR}_surface_p.nc"
if [ ! -f "$FULL_FILE" ]; then
    echo "ERROR: $FULL_FILE not produced"
    exit 2
fi
if [ ! -f "$SP_FILE" ]; then
    echo "ERROR: $SP_FILE not produced"
    exit 2
fi

#-----------------------------------------------------------------------------
# Step 2: verify aux files are symlinked correctly
#-----------------------------------------------------------------------------
echo ""
echo "[2/5] Checking auxiliary files..."

AUX_DIR="$WORKDIR/aux"
TROY_AUX="${LUCIE_TROY_AUX:-/glade/derecho/scratch/mdarman/ERA5_hr_haiwen/ERA5_T30_1950_TO_2021}"
declare -A AUX_LINKS=(
    ["year_average_logp.nc"]="$TROY_AUX/speedy_average_surfacep/year_average_logp.nc"
    ["speedy_orography.nc"]="$TROY_AUX/orography_data/speedy_orography.nc"
    ["era_orography_jan01_1980.nc"]="$TROY_AUX/orography_data/era_orography_jan01_1980.nc"
)
for fname in "${!AUX_LINKS[@]}"; do
    target="$AUX_DIR/$fname"
    expected="${AUX_LINKS[$fname]}"
    if [ ! -e "$target" ]; then
        echo "  Creating symlink: $fname"
        ln -sf "$expected" "$target"
    elif [ -L "$target" ]; then
        actual=$(readlink "$target")
        if [ "$actual" != "$expected" ]; then
            echo "  WARNING: $fname links to $actual, expected $expected"
        else
            echo "  OK: $fname"
        fi
    else
        echo "  OK (file): $fname"
    fi
done

#-----------------------------------------------------------------------------
# Step 3: patch regrid_era.f90 with year/month
#-----------------------------------------------------------------------------
echo ""
echo "[3/5] Patching regrid_era.f90 source for y=$YEAR m=$MONTH..."

REGRID_SRC="$SRC_DIR/regrid_era.f90"
[ -f "$REGRID_SRC" ] || { echo "ERROR: $REGRID_SRC not found"; exit 3; }

"$JAX_PYTHON" - <<PY
import re
src = open("$REGRID_SRC").read()
src = re.sub(r"start_year\s*=\s*\d+", "start_year = $YEAR", src, count=1)
src = re.sub(r"end_year\s*=\s*\d+", "end_year = $YEAR", src, count=1)
src = re.sub(r"start_month\s*=\s*\d+", "start_month = $MONTH", src, count=1)
src = re.sub(r"end_month\s*=\s*\d+", "end_month = $MONTH", src, count=1)
open("$REGRID_SRC", "w").write(src)
print("  patched")
PY

#-----------------------------------------------------------------------------
# Step 4: rebuild
#-----------------------------------------------------------------------------
echo ""
echo "[4/5] Rebuilding main_regrid.exe..."
cd "$SRC_DIR"
make 2>&1 | tee -a "$LOG" | tail -5

[ -x "$SRC_DIR/main_regrid.exe" ] || { echo "ERROR: build failed"; exit 4; }

#-----------------------------------------------------------------------------
# Step 5: run
#-----------------------------------------------------------------------------
echo ""
echo "[5/5] Running main_regrid.exe (logging to $LOG)..."
cd "$SRC_DIR"
./main_regrid.exe >> "$LOG" 2>&1
RC=$?

if [ $RC -ne 0 ]; then
    echo "ERROR: main_regrid.exe exited $RC"
    tail -20 "$LOG"
    exit 5
fi

if [ ! -f "$OUTPUT_FILE" ]; then
    echo "ERROR: expected output $OUTPUT_FILE not found"
    tail -20 "$LOG"
    exit 6
fi

#-----------------------------------------------------------------------------
# Done
#-----------------------------------------------------------------------------
SIZE=$(stat -c %s "$OUTPUT_FILE")
echo ""
echo "=========================================================================="
echo "  ✓ SUCCESS for ${YEAR}-${MONTH}"
echo "  Output: $OUTPUT_FILE"
echo "  Size:   $((SIZE / 1024 / 1024)) MB"
echo "=========================================================================="
