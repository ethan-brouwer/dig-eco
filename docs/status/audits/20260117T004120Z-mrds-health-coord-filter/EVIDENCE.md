# Evidence: 20260117T004120Z-mrds-health-coord-filter

## Scope
Filter MRDS rows with missing or blank coordinates to prevent
`ee.Number(null)` errors when building point geometries.

## Pre-state
- git status recorded in evidence/git-status-pre.txt

## Changes
- Updated `gee/groundwork/MDRS_health.js`
  - Filtered out null/blank coordinate rows
  - Parsed coordinate strings before geometry creation

## Post-state
- git status recorded in evidence/git-status-post.txt
- git diff --stat recorded in evidence/git-diff-stat-post.txt

## Validation
No runtime tests executed (GEE script only).
