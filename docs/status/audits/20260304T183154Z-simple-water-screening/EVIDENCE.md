# Evidence: 20260304T183154Z-simple-water-screening

## Scope
Add a GitHub-ready GEE script for simple upstream/downstream water screening at
San Sebastian Mine using recent Sentinel-2 imagery and conservative water-only
proxies.

## Pre-state
- git status recorded in `evidence/git-status-pre.txt`

## Changes
- Added `gee/groundwork/san_sebastian_upstream_downstream_water_simple.js`
  with:
  - recent composite only
  - JRC + optical water mask
  - one upstream control and three downstream bins
  - `TSS_PROXY`, `NDTI`, and `RED_GREEN` reach summaries
  - QA layers and quick charts for visual verification in GEE

## Post-state
- git status recorded in `evidence/git-status-post.txt`
- git diff --stat recorded in `evidence/git-diff-stat-post.txt`

## Validation
Script-only step. No runtime tests executed locally.
