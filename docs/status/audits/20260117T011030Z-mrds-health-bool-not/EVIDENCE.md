# Evidence: 20260117T011030Z-mrds-health-bool-not

## Scope
Replace `.not()` on `ee.Algorithms.IsEqual` with `ee.Algorithms.Not`
to avoid boolean method errors.

## Pre-state
- git status recorded in evidence/git-status-pre.txt

## Changes
- Updated `gee/groundwork/MDRS_health.js`
  - Corrected boolean negation for coordinate validation

## Post-state
- git status recorded in evidence/git-status-post.txt
- git diff --stat recorded in evidence/git-diff-stat-post.txt

## Validation
No runtime tests executed (GEE script only).
