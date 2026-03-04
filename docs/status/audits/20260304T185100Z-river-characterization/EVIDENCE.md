# Evidence: 20260304T185100Z-river-characterization

## Scope
Add a dedicated San Sebastian river characterization script focused only on:
- river mask troubleshooting
- water-only proxy visualization
- transparent QA layers for why the river may disappear

## Pre-state
- git status recorded in `evidence/git-status-pre.txt`

## Changes
- Added `Phase I/gee_scripts/turbidity/san_sebastian_river_characterization.js`
- Added explicit debug layers for:
  - JRC water occurrence
  - JRC river seed
  - JRC river corridor
  - optical water mask
  - final river analysis mask
  - disagreement masks between JRC and optical water

## Post-state
- git status recorded in `evidence/git-status-post.txt`
- git diff --stat recorded in `evidence/git-diff-stat-post.txt`

## Validation
Script-only step. No local runtime test executed. Final validation must occur
in the GEE Code Editor by checking whether the final river analysis mask aligns
with the visible river.
