# Evidence: 20260116T234430Z-mrds-health-empty-months

## Scope
Skip empty monthly composites to avoid no-band errors in MRDS health workflow.

## Pre-state
- git status recorded in evidence/git-status-pre.txt

## Changes
- Updated `gee/groundwork/MDRS_health.js` to drop empty months.

## Post-state
- git status recorded in evidence/git-status-post.txt
- git diff --stat recorded in evidence/git-diff-stat-post.txt

## Validation
Script-only step. No runtime tests executed.
