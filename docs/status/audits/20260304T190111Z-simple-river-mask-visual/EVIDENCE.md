# Evidence: 20260304T190111Z-simple-river-mask-visual

## Scope
Add a minimal Phase I river-mask visualization script focused only on quick
visual QA of channel coverage.

## Pre-state
- git status recorded in `evidence/git-status-pre.txt`

## Changes
- Added `Phase I/gee_scripts/turbidity/san_sebastian_river_mask_visual_simple.js`
  with:
  - true/false-color context
  - optical water occurrence fraction
  - optical instant mask
  - final combined river mask
  - minimal threshold controls

## Post-state
- git status recorded in `evidence/git-status-post.txt`
- git diff --stat recorded in `evidence/git-diff-stat-post.txt`

## Validation
Script-only step. No local runtime test executed. Final validation should be
done in GEE by visual comparison against the true-color/false-color channel.
