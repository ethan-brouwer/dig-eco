# Evidence: 20260304T184122Z-river-visualization

## Scope
Patch the simple San Sebastian characterization script to:
- show the river analysis area explicitly
- remove reliance on inferred upstream/downstream bins for this step
- add a water-only proxy composite for visual inspection along the channel

## Pre-state
- git status recorded in `evidence/git-status-pre.txt`

## Changes
- Updated `gee/groundwork/san_sebastian_site_characterization_simple.js` to add:
  - `Coarse river context (HydroSHEDS)`
  - `Current water mask (optical + JRC)`
  - `River analysis mask`
  - `QA excluded from river analysis mask`
  - `Water proxy composite (TSS / NDTI / red-green)`

## Post-state
- git status recorded in `evidence/git-status-post.txt`
- git diff --stat recorded in `evidence/git-diff-stat-post.txt`

## Validation
Script-only step. No local runtime test executed. Final validation must occur in
the GEE Code Editor by checking whether the river analysis mask aligns with the
visibly interpretable channel.
