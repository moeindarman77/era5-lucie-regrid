# Vertical Interpolation: ERA5 Pressure Levels to SPEEDY Sigma Levels

## Overview

After horizontal regridding to the T30 Gaussian grid, the data must be
vertically interpolated from ERA5's 28 pressure levels to SPEEDY's 8 sigma
levels. This uses 1D B-spline interpolation with climatological surface
pressure.

## Input

28 ERA5 pressure levels (in hPa), a subset of ERA5's 37 native levels:

```
20, 30, 50, 70, 100, 125, 150, 175, 200, 225, 250, 300, 350, 400, 450,
500, 550, 600, 650, 700, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975
```

Note: levels 1, 2, 3, 5, 7, 10, 15, and 1000 hPa are excluded to match
Troy's original specification.

## Output

8 SPEEDY sigma levels:

| Index | Sigma (σ) | Approximate altitude |
|-------|-----------|---------------------|
| 0     | 0.025     | ~35 km (stratosphere) |
| 1     | 0.095     | ~17 km |
| 2     | 0.200     | ~12 km |
| 3     | 0.340     | ~8 km |
| 4     | 0.510     | ~5.5 km |
| 5     | 0.685     | ~3 km |
| 6     | 0.835     | ~1.5 km |
| 7     | 0.950     | ~0.5 km (near surface) |

Sigma coordinates are defined as σ = P / P_surface, where P is the
pressure at the level and P_surface is the surface pressure.

## Method: Climatological Sigma

The vertical interpolation is performed by `era_p_level_to_speedy_p_level`
in `fortran/mod_interp2d.f90` (line ~305).

For each horizontal grid point (i, j):

1. **Compute target pressures**: P_target(k) = σ(k) × P_surface_clim(i, j)
   where P_surface_clim is the **annual-mean** surface pressure from
   `year_average_logp.nc`
2. **Interpolate**: Use 1D B-spline (`bspline_interface1d`) to interpolate
   the ERA5 profile at the 28 pressure levels to the 8 target pressures

This is "climatological sigma" — the sigma-to-pressure conversion uses a
fixed annual-mean surface pressure field rather than the instantaneous surface
pressure at each timestep. This matches Troy's methodology and produces
bit-perfect reproduction of the training data.

### Why Climatological (Not True) Sigma?

The pipeline supports both:
- `era_p_level_to_speedy_p_level` — uses annual-mean SP (**this is what's used**)
- `era_p_level_to_speedy_sigma_level` — uses instantaneous SP (commented out)

The climatological approach was validated as bit-perfect against Troy's
training data. Using true sigma would give physically more accurate
interpolation but would NOT reproduce the training data.

## The `year_average_logp.nc` File

This auxiliary file contains the annual-mean log(surface pressure) on the
T30 Gaussian grid (96×48). It was computed by Troy from the full 1950-2021
ERA5 record.

The pipeline:
1. Reads `average_logp` from the file (shape: 96×48×1)
2. Exponentiates: `P_surface = exp(average_logp)`
3. Converts to Pa: `P_surface = P_surface × 1000`

**Critical**: This file MUST be Troy's actual file from
`ERA5_T30_1950_TO_2021/speedy_average_surfacep/year_average_logp.nc`.
Reconstructing it from LUCIE h5 `logp` values produces a slightly different
distribution that causes ~1 K temperature errors after interpolation.

## B-spline Order for Vertical Interpolation

The 1D vertical interpolation uses `kx = 3` (in bspline-fortran terms).

Using the same FITPACK convention as horizontal interpolation:

| FITPACK convention | bspline-fortran `kx` | Actual spline degree |
|--------------------|----------------------|---------------------|
| k=3 (cubic)        | kx=3                 | degree 2 (piecewise quadratic) |

Troy's paper describes "1D cubic B-spline" for the vertical interpolation.
The `kx=3` parameter produces bit-perfect reproduction.

## Key Source Files

- `fortran/mod_interp2d.f90`:
  - `speedy_sigma_levels` parameter (line ~8): the 8 sigma values
  - `era_p_level_to_speedy_p_level` (line ~305): climatological sigma interpolation
  - `bspline_interface1d` (line ~145): 1D B-spline wrapper
- `fortran/regrid_era.f90`: calls the vertical interpolation at line ~211
- `aux/year_average_logp.nc`: annual-mean surface pressure field
