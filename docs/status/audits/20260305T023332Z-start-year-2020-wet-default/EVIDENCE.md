# Evidence: 20260305T023332Z-start-year-2020-wet-default

## Scope
Set active turbidity comparison workflows to analyze from year 2020 through
present and keep wet-season filtering enabled by default.

## Pre-state
- git status recorded in `evidence/git-status-pre.txt`

## Changes
- Updated `Phase I/gee_scripts/turbidity/san_sebastian_ssrivers_monthly_simple.js`
- Updated `Phase I/gee_scripts/turbidity/san_sebastian_updown_lines_compare_simple.js`
- Replaced rolling-year window with fixed `analysisStartYear = 2020`
- Added explicit print of start year in script output

## Post-state
- git status recorded in `evidence/git-status-post.txt`
- git diff --stat recorded in `evidence/git-diff-stat-post.txt`

## Validation
Script-only step. No local runtime tests executed. Final validation is in GEE
by confirming analysis window and charts include years from 2020 onward.
