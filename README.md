# flat-entropy-envelope

Replace the envelope of a MESA stellar profile with a **constant-entropy (isentropic) envelope** in hydrostatic equilibrium, for use as initial conditions in [Phantom](https://github.com/danieljprice/phantom) SPH simulations.

The material inside a chosen core radius is kept from the MESA profile. Outside it, the code integrates the equations of stellar structure outward on the MESA mass grid at fixed entropy, and adjusts that entropy until the surface pressure matches a target value. The result is written as a MESA-format profile that Phantom's star setup can read.

Typical use: build a convective envelope with no superadiabatic or radiative surface layer, so that the stellar envelope stays put when relaxed and evolved in SPH.

## Contents

| File | Purpose |
|---|---|
| `fixed_S_envelope.f90` | Module `fixedSprofile`: builds the isentropic envelope (`calc_rho_and_pres`, `one_shot`). |
| `write_fixedSprofile.f90` | Driver program: reads the MESA profile, calls the envelope solver, computes T and u, writes the output profile. |
| `write_uTprofile.f90` | Utility: recomputes temperature, internal energy and (optionally) sound speed of a profile for a different equation of state (e.g. ideal gas, `ieos = 2`). |
| `compare_profiles.py` | Jupyter-cell script (VS Code `# %%` cells) that plots the new profile against the MESA profile. |
| `Makefile` | Copies the sources into the Phantom tree and builds them there. |
| `Makefile_phantom` | Phantom build makefile with extra targets for the programs above. |
| `phantom/` | Phantom, as a git submodule. |
| `profile40.data` | Example MESA profile (12 M☉ red supergiant, R ≈ 619 R☉). |

## Requirements

- A Fortran compiler supported by Phantom (tested with `gfortran`).
- Phantom, checked out as the `phantom` submodule:
  ```
  git submodule update --init
  ```
- For plotting: Python 3 with `numpy` and `matplotlib`, and VS Code with the Jupyter extension (or any Jupyter environment).

## Building

From the repository root:

```
make fixedSenvelope
```

This copies `Makefile_phantom`, `fixed_S_envelope.f90` and `write_fixedSprofile.f90` into the Phantom tree, builds the `fixedSprofile` target there against Phantom's libraries (with `SETUP=star`), and copies the `fixedSprofile` executable back to the repository root.

The other targets work the same way:

```
make writeuT      # builds write_uTprofile
make clean        # cleans the Phantom build directory
```

Because the sources are copied into `phantom/src/setup/`, always rebuild through this `Makefile` after editing them. Do not edit the copies inside `phantom/`.

## Running

```
./fixedSprofile
```

The program:

1. Reads `./profile40.data`.
2. Prints the core radius, target surface pressure, the composition (X, Z, μ) taken from the middle of the envelope, and the initial entropy guess, then **waits for Enter**.
3. Shoots on the envelope entropy, printing one line per shot (shot number, `ok`/`failed`, entropy, surface pressure, radius reached).
4. Prints the obtained surface pressure and radius.
5. Computes temperature and specific internal energy everywhere from ρ and P.
6. Writes `fixedSprofile.dat`.

### Settings

The input and output file names are set in `write_fixedSprofile.f90`. Everything else is in the `USER SETTINGS` block of `calc_rho_and_pres` in `fixed_S_envelope.f90`:

| Variable | Meaning | Example |
|---|---|---|
| `rcore` | Core radius. Material inside it is kept from MESA. | `18.5 * solarr` |
| `surfpres` | Target surface pressure (dyn cm⁻²). | `600.` |
| `Sfac` | Entropy step (k_B per particle) used while bracketing the solution. | `0.05` |
| `tol` | Relative tolerance on the surface pressure. | `1d-3` |
| `maxiter` | Maximum number of shots. | `200` |
| `ieos` | Phantom equation of state used for the envelope. | `12` (ideal gas + radiation) |
| `ientropy` | Entropy definition passed to Phantom's `entropy()`. | `2` (gas + radiation) |

Recompile after changing any of them.

## Method

### Equations

Outside `rcore` the envelope satisfies hydrostatic equilibrium and mass conservation in Lagrangian form,

$$\frac{dP}{dm} = -\frac{G m}{4\pi r^4}, \qquad \frac{dr}{dm} = \frac{1}{4\pi r^2 \rho},$$

closed by the condition that the specific entropy is the same everywhere in the envelope, $s(\rho, P) = S_c$. Entropy is evaluated with Phantom's `entropy()` (ideal gas plus radiation for `ientropy = 2`), with a single mean molecular weight μ taken from the MESA profile halfway between `rcore` and the surface.

### Grid

The envelope is integrated on the MESA mass grid, so the mass coordinate and total mass are exactly those of the input profile. The grid follows the MESA staggered convention, numbered `j = 0 … N` from the core boundary to the surface:

```
     cell 0   |  cell 1   |  cell 2   | ... |  cell N  |
  (in core) face 0      face 1      face 2  ...       face N = surface
             r(0) ≈ rcore
```

| Array | Bounds | Definition |
|---|---|---|
| `m, r` | `0:N` | mass and radius at face `j` (outer edge of cell `j`) |
| `rho, pres` | `0:N` | averages over cell `j`, which lies between faces `j-1` and `j` |
| `dm` | `1:N` | `dm(j) = m(j) − m(j−1)`, mass of cell `j` |
| `dm_c` | `0:N−1` | `dm_c(j) = ½(m(j+1) − m(j−1))`, mass between the centres of cells `j` and `j+1` |

Point `j = 0` is the MESA point nearest `rcore`. Its `r`, `rho` and `pres` are taken from MESA and form the inner boundary condition. MESA point `core_idx + j` maps to envelope point `j`.

### Integration (`one_shot`)

For a given entropy `S_c`, each step moves outward across face `j`:

1. Pressure across the face: `pres(j+1) = pres(j) − G m(j) dm_c(j) / (4π r(j)⁴)`
2. Density from the equation of state: `rho(j+1) = rho(pres(j+1), S_c)`
3. Radius of the next face: `r(j+1)³ = r(j)³ + 3 dm(j+1) / (4π rho(j+1))`

A shot fails if the pressure or density becomes negative before the surface is reached.

### Shooting on entropy

The surface pressure `pres(N)` increases with `S_c`. A failed shot means the entropy is too low: the envelope is too dense and runs out of pressure before all the mass is placed. The solver:

1. Starts from the entropy of the MESA profile at mid-envelope.
2. Steps `S_c` up or down by `Sfac` until it has one entropy that gives too low a surface pressure (or fails) and one that gives too high.
3. Bisects between them until `|pres(N) − surfpres| < tol · surfpres`.

It stops with an error if the bracket collapses without converging or `maxiter` is reached.

### Output

`fixedSprofile.dat` is a MESA-format profile in cgs units, ordered from surface to centre, with columns

```
Mass  Pressure  Temperature  Radius  Density  Eint
```

Inside `rcore` it contains the original MESA ρ, P and r. Temperature and specific internal energy are recomputed everywhere from ρ and P with the chosen `ieos`. Composition is not written, so give Phantom the same X, Z and μ that the code printed when it started.


## Limitations

- **The stellar radius is not controlled.** Only the surface pressure is matched. With the core boundary fixed to the MESA values, an isentropic envelope holding the same mass is generally more extended than the MESA star, because the real outer envelope has lower entropy and higher μ.
- **Single μ.** The envelope uses one mean molecular weight. Where hydrogen recombines near the surface, the real μ is higher.
- **Inner boundary.** `rho(0)` and `pres(0)` come from MESA and are not on the adiabat `S_c`, so there is a small density discontinuity at the core boundary.
- **The core radius snaps to the MESA grid.** The boundary is the MESA point nearest `rcore`. Avoid placing `rcore` inside a composition discontinuity (in `profile40.data`, the edge of the hydrogen-burning shell is at 18.37 R☉).