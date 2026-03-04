# Evidence: 20260304T191323Z-wet-season-river-mask

## Scope
Update the simple river-mask visual QA script to support wet-season-only image
selection so river visibility is less affected by dry-month channel contraction.

## Pre-state
- git status recorded in `evidence/git-status-pre.txt`

## Changes
- Updated `Phase I/gee_scripts/turbidity/san_sebastian_river_mask_visual_simple.js`
  to:
  - add `useWetSeasonFilter`
  - add `wetSeasonMonths`
  - add month metadata and month filtering to Sentinel-2 collection build
  - print active wet-season filter settings

## Post-state
- git status recorded in `evidence/git-status-post.txt`
- git diff --stat recorded in `evidence/git-diff-stat-post.txt`

## Validation
Script-only step. No local runtime test executed. Final validation should be
done in GEE by verifying that wet-season filtering increases plausible river
coverage while avoiding major non-river false positives.
