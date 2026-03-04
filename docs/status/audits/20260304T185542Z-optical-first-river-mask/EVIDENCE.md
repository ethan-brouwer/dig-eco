# Evidence: 20260304T185542Z-optical-first-river-mask

## Scope
Update the Phase I river characterization script so the river analysis mask is
optical-first, with HydroSHEDS and JRC retained as comparison and QA layers
rather than hard gating filters.

## Pre-state
- git status recorded in `evidence/git-status-pre.txt`

## Changes
- Updated `Phase I/gee_scripts/turbidity/san_sebastian_river_characterization.js`
  to:
  - make `opticalWater` the primary `riverAnalysisMask`
  - add `HydroSHEDS river seed` and `HydroSHEDS river corridor`
  - add `Hydro-constrained river mask` and `JRC-constrained river mask`
  - add disagreement QA layers for optical vs HydroSHEDS and optical vs JRC

## Post-state
- git status recorded in `evidence/git-status-post.txt`
- git diff --stat recorded in `evidence/git-diff-stat-post.txt`

## Validation
Script-only step. No local runtime test executed. Final validation must occur in
the GEE Code Editor by checking whether the optical-first river mask tracks the
visible channel better than the constrained alternatives.
