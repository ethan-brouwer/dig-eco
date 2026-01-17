# Evidence: 20260117T010230Z-mrds-health-monthly-xaxis

## Scope
Restore monthly dry-season series and harden coordinate parsing while
keeping existing map/selector behavior intact.

## Pre-state
- git status recorded in evidence/git-status-pre.txt

## Changes
- Updated `gee/groundwork/MDRS_health.js`
  - Monthly dry-season composites for charting
  - Readable x-axis and chart layout
  - Regex-based coordinate validation and filtering

## Post-state
- git status recorded in evidence/git-status-post.txt
- git diff --stat recorded in evidence/git-diff-stat-post.txt

## Validation
No runtime tests executed (GEE script only).
