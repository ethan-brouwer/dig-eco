# Evidence: 20260117T012300Z-mrds-health-no-not

## Scope
Remove `.not()` usage to avoid reliance on `ee.Algorithms.Not` in masks.

## Pre-state
- git status recorded in evidence/git-status-pre.txt

## Changes
- Updated `gee/groundwork/MDRS_health.js`
  - Replace image `.not()` calls with `.eq(0)` masks

## Post-state
- git status recorded in evidence/git-status-post.txt
- git diff --stat recorded in evidence/git-diff-stat-post.txt

## Validation
No runtime tests executed (GEE script only).
