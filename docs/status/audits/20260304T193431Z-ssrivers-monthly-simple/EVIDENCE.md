# Evidence: 20260304T193431Z-ssrivers-monthly-simple

## Scope
Add a simplified SSrivers-based workflow that:
- buffers user-drawn river lines
- maps water-only index proxies inside the corridor
- outputs monthly charts for key proxies

## Pre-state
- git status recorded in `evidence/git-status-pre.txt`

## Changes
- Added `Phase I/gee_scripts/turbidity/san_sebastian_ssrivers_monthly_simple.js`
  with:
  - manual `SSrivers` geometry handling
  - buffered corridor mask
  - river mask from NDWI/MNDWI inside corridor
  - water-only proxy layers (`TSS_PROXY`, `NDTI`, `RED_GREEN`)
  - monthly chart outputs for each proxy

## Post-state
- git status recorded in `evidence/git-status-post.txt`
- git diff --stat recorded in `evidence/git-diff-stat-post.txt`

## Validation
Script-only step. No local runtime test executed. Final validation occurs in GEE
by visually checking corridor/mask alignment and reviewing monthly charts.
