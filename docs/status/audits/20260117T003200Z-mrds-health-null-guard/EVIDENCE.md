# Evidence: 20260117T003200Z-mrds-health-null-guard

## Scope
Guard against null values during MRDS index normalization to prevent
`ee.Number(null)` errors.

## Pre-state
- git status recorded in evidence/git-status-pre.txt

## Changes
- Updated `gee/groundwork/MDRS_health.js`
  - Added safe numeric fallback handling
  - Hardened max and per-feature normalization

## Post-state
- git status recorded in evidence/git-status-post.txt
- git diff --stat recorded in evidence/git-diff-stat-post.txt

## Validation
No runtime tests executed (GEE script only).
