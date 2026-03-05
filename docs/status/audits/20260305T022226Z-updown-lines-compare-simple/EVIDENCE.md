# Evidence: 20260305T022226Z-updown-lines-compare-simple

## Scope
Add a simple upstream/downstream comparison workflow using user-drawn line
imports (`upstream1`, `downstream1`) and `impact_point`.

## Pre-state
- git status recorded in `evidence/git-status-pre.txt`

## Changes
- Added `Phase I/gee_scripts/turbidity/san_sebastian_updown_lines_compare_simple.js`
  with:
  - buffered corridors for upstream/downstream lines
  - monthly proxy summaries by reach (`TSS_PROXY`, `NDTI`, `RED_GREEN`)
  - side-by-side monthly charts by reach type
  - monthly downstream-minus-upstream difference charts

## Post-state
- git status recorded in `evidence/git-status-post.txt`
- git diff --stat recorded in `evidence/git-diff-stat-post.txt`

## Validation
Script-only step. No local runtime tests executed. Final validation should be
done in GEE by checking line/corridor alignment and chart behavior.
