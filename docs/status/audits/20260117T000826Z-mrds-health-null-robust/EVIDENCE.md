# Evidence: 20260117T000826Z-mrds-health-null-robust

## Scope
Guard null values in MRDS health normalization with explicit null checks.

## Pre-state
- git status recorded in evidence/git-status-pre.txt

## Changes
- Updated `gee/groundwork/MDRS_health.js` to treat nulls as 0 via IsEqual checks.

## Post-state
- git status recorded in evidence/git-status-post.txt
- git diff --stat recorded in evidence/git-diff-stat-post.txt

## Validation
Script-only step. No runtime tests executed.
