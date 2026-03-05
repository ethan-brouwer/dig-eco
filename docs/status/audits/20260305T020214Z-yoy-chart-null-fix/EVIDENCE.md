# Evidence: 20260305T020214Z-yoy-chart-null-fix

## Scope
Fix year-over-year chart generation error caused by setting `null` values in
dictionary pivot logic.

## Pre-state
- git status recorded in `evidence/git-status-pre.txt`

## Changes
- Updated `Phase I/gee_scripts/turbidity/san_sebastian_ssrivers_monthly_simple.js`
  to:
  - replace dictionary-pivot overlay builder with long-form grouped charts
    (`ui.Chart.feature.groups`)
  - use `*_chart` fields and skip no-data points naturally
  - keep Jan-Dec month labeling via numeric month ticks

## Post-state
- git status recorded in `evidence/git-status-post.txt`
- git diff --stat recorded in `evidence/git-diff-stat-post.txt`

## Validation
Script-only step. No local runtime tests executed. Final validation in GEE
should confirm YOY charts render without dictionary-null errors.
