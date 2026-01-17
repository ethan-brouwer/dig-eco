# Evidence: 20260116T195500Z-mrds-health-qol

## Scope
Fix normalization mapping and improve map visibility for MRDS health workflow.

## Pre-state
- git status recorded in evidence/git-status-pre.txt

## Changes
- Updated `gee/groundwork/MDRS_health.js` to use ee.List for normalization.
- Ensured MRDS sites and El Salvador outline are visible by default.

## Post-state
- git status recorded in evidence/git-status-post.txt
- git diff --stat recorded in evidence/git-diff-stat-post.txt

## Validation
Script-only step. No runtime tests executed.
