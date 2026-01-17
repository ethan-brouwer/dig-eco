# Evidence: 20260116T194433Z-mrds-health-fix

## Scope
Update MRDS health workflow to use monthly dry-season medians and robust normalization.

## Pre-state
- git status recorded in evidence/git-status-pre.txt

## Changes
- Updated `gee/groundwork/MDRS_health.js` for per-month composites and safe normalization.

## Post-state
- git status recorded in evidence/git-status-post.txt
- git diff --stat recorded in evidence/git-diff-stat-post.txt

## Validation
Script-only step. No runtime tests executed.
