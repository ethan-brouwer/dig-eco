# Evidence: 20260304T185542Z-optical-occurrence-river-mask

## Scope
Improve the Phase I river characterization mask to avoid near-empty results on
narrow/intermittent channels by using optical water occurrence across the full
image collection (not only a median composite threshold).

## Pre-state
- git status recorded in `evidence/git-status-pre.txt`

## Changes
- Updated `Phase I/gee_scripts/turbidity/san_sebastian_river_characterization.js`
  to:
  - add optical water occurrence fraction (`OPTICAL_WATER_FRACTION`)
  - create an occurrence mask (`OPTICAL_WATER_OCCURRENCE`)
  - define primary optical mask as occurrence OR instant median mask
  - add hydro-assisted comparison mask (`HYDRO_ASSISTED_RIVER_MASK`)
  - add additional diagnostics and visualization layers for occurrence behavior

## Post-state
- git status recorded in `evidence/git-status-post.txt`
- git diff --stat recorded in `evidence/git-diff-stat-post.txt`

## Validation
Script-only step. No local runtime test executed. Final validation must occur
in GEE by checking that optical occurrence recovers plausible river pixels
before applying HydroSHEDS/JRC constraints.
