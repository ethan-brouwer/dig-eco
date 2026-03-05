# Evidence: 20260305T014919Z-corridor-direct-monthly

## Scope
Update the simplified SSrivers monthly workflow to support corridor-direct
monthly calculations (50 m buffer) for more stable month-to-month output when
water-pixel masks are sparse.

## Pre-state
- git status recorded in `evidence/git-status-pre.txt`

## Changes
- Updated `Phase I/gee_scripts/turbidity/san_sebastian_ssrivers_monthly_simple.js`
  to:
  - default to `corridorBufferMeters = 50`
  - add `useCorridorDirectly = true`
  - use corridor+observation area as monthly analysis mask when enabled
  - keep monthly QA flags (`corridor_direct`, `monthly_mask`, `fallback_reference_mask`, `no_images`)
  - expose additional mask layers for QA

## Post-state
- git status recorded in `evidence/git-status-post.txt`
- git diff --stat recorded in `evidence/git-diff-stat-post.txt`

## Validation
Script-only step. No local runtime tests executed. Final validation is in GEE
by checking monthly chart continuity and QA flags.
