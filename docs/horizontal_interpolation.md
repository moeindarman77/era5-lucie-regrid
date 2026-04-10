# Horizontal Interpolation: ERA5 to T30 Gaussian Grid

## Overview

The horizontal regridding takes ERA5 data on a regular latitude-longitude grid
and interpolates it onto the SPEEDY T30 Gaussian grid using 2D B-spline
interpolation.

## Grids

**Input**: ERA5 at 3.0° resolution
- 120 longitudes: 0°, 3°, 6°, ..., 357°
- 61 latitudes: 90°N, 87°N, ..., 87°S, 90°S

**Output**: SPEEDY T30 Gaussian grid
- 96 longitudes: 0°, 3.75°, 7.5°, ..., 356.25°
- 48 Gaussian latitudes (S→N in output)

## Method

The interpolation uses the `bspline_interface` subroutine in
`fortran/mod_interp2d.f90`, which wraps the
[bspline-fortran](https://github.com/jacobwilliams/bspline-fortran) library's
`db2ink` (setup) and `db2val` (evaluation) routines.

For each 2D field at each timestep:

1. **Set up** the B-spline representation via `db2ink(oldx, nx, oldy, ny, var, kx, ky, iknot, tx, ty, bcoef, iflag)`
2. **Evaluate** at each target grid point via `db2val(newx(i), newy(j), ..., bcoef, val, ...)`

## B-spline Order (the FITPACK off-by-one)

The code uses `kx = ky = 2`.

This is a naming convention subtlety between FITPACK and bspline-fortran:

| FITPACK convention | bspline-fortran `kx` | Actual spline degree |
|--------------------|----------------------|---------------------|
| k=1 (linear)       | kx=1                 | degree 0 (piecewise constant) |
| k=2 (quadratic)    | kx=2                 | degree 1 (piecewise linear) |
| k=3 (cubic)        | kx=3                 | degree 2 (piecewise quadratic) |

Troy's paper describes "quadratic B-spline" (FITPACK k=2). The original code
inherited this convention, and `kx=ky=2` is what produces bit-perfect
reproduction of the LUCIE training data.

## Latitude Handling

ERA5 data arrives with latitudes ordered North→South (90°N to 90°S). The
regridder needs South→North ordering internally:

1. ERA5 latitudes are **flipped** in `regrid_era.f90` (line ~150) before being
   passed to the interpolation routines
2. The output is written with a `gridy:1:-1` index reversal (line ~206),
   which flips the result back — but because both the input flip and output
   flip cancel, the final output ends up in **S→N order**, matching the LUCIE
   training data convention

## Variables Interpolated

Each of the 4 prognostic variables (T, U, V, Q) at all 28 pressure levels,
plus surface pressure (SP), are independently interpolated from the ERA5
grid to the T30 grid using this 2D B-spline method.

The horizontal interpolation happens *before* the vertical interpolation to
sigma levels. The vertically interpolated result is then the final output.

## Extrapolation

Extrapolation is enabled (`extrap = .True.`) for grid points that fall
slightly outside the ERA5 domain. In practice, this mainly affects the
longitude wrapping (ERA5 goes to 357° while SPEEDY goes to 356.25°).

## Key Source Files

- `fortran/mod_interp2d.f90`: `bspline_interface` subroutine (line ~198)
- `fortran/regrid_era.f90`: main loop calling the interpolation (line ~200+)
- `fortran/mod_utilities.f90`: grid constants (`gridx=96`, `gridy=48`, `speedy_lat`)
