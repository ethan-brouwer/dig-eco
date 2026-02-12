# Evidence: 20260117T011700Z-mrds-health-or-if

## Scope
Replace boolean negation with `ee.Algorithms.Or` + `ee.Algorithms.If`
to avoid unsupported `.Not` usage.

## Pre-state
- git status recorded in evidence/git-status-pre.txt

## Changes
- Updated `gee/groundwork/MDRS_health.js`
  - Compute valid_coord with Or/If instead of Not

## Post-state
- git status recorded in evidence/git-status-post.txt
- git diff --stat recorded in evidence/git-diff-stat-post.txt

## Validation
No runtime tests executed (GEE script only).
