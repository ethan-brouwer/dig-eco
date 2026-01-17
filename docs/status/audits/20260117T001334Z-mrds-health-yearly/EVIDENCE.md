# Evidence: 20260117T001334Z-mrds-health-yearly

## Scope
Adjust MRDS health charting to plot a single dry-season point per year
and normalize each index by its maximum over the analysis period.

## Pre-state
- git status recorded in evidence/git-status-pre.txt

## Changes
- Updated `gee/groundwork/MDRS_health.js`
  - Yearly dry-season mean per site
  - Max-over-period normalization
  - Cleaner chart options

## Post-state
- git status recorded in evidence/git-status-post.txt
- git diff --stat recorded in evidence/git-diff-stat-post.txt

## Validation
No runtime tests executed (GEE script only).
