# Evidence: 20260117T000155Z-mrds-health-null-guard

## Scope
Guard null values in MRDS health normalization.

## Pre-state
- git status recorded in evidence/git-status-pre.txt

## Changes
- Updated `gee/groundwork/MDRS_health.js` to default null values to 0 before normalization.

## Post-state
- git status recorded in evidence/git-status-post.txt
- git diff --stat recorded in evidence/git-diff-stat-post.txt

## Validation
Script-only step. No runtime tests executed.
