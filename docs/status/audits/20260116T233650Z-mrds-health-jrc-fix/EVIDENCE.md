# Evidence: 20260116T233650Z-mrds-health-jrc-fix

## Scope
Fix JRC Global Surface Water asset reference in MRDS health workflow.

## Pre-state
- git status recorded in evidence/git-status-pre.txt

## Changes
- Updated `gee/groundwork/MDRS_health.js` to use GlobalSurfaceWater/occurrence.

## Post-state
- git status recorded in evidence/git-status-post.txt
- git diff --stat recorded in evidence/git-diff-stat-post.txt

## Validation
Script-only step. No runtime tests executed.
