# Evidence: 20260304T192641Z-ssrivers-buffer-indices

## Scope
Add a script that uses user-drawn `SSrivers` geometry to define a buffered
river corridor and run water-index screening only inside that corridor.

## Pre-state
- git status recorded in `evidence/git-status-pre.txt`

## Changes
- Added `Phase I/gee_scripts/turbidity/san_sebastian_ssrivers_buffer_indices.js`
  with:
  - support for `SSrivers` imported as Geometry/Feature/FeatureCollection
  - buffered corridor mask from manual river geometry
  - wet-season Sentinel-2 option
  - corridor-constrained optical river mask
  - corridor-constrained water proxy layers: `TSS_PROXY`, `NDTI`, `RED_GREEN`

## Post-state
- git status recorded in `evidence/git-status-post.txt`
- git diff --stat recorded in `evidence/git-diff-stat-post.txt`

## Validation
Script-only step. No local runtime test executed. Final validation in GEE
should confirm corridor geometry alignment and river-mask coverage.
