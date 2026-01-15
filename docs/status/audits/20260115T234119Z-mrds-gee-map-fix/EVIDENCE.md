# Evidence: 20260115T234119Z-mrds-gee-map-fix

## Scope
Fix MRDS GEE script for CSV imports (remove 3-arg get) and allow
configurable field names for labels and commodity styling.

## Pre-state
- git status recorded in evidence/git-status-pre.txt

## Changes
- Updated `mrds_gee_map.js` for safe property access and CSV field config.

## Post-state
- git status recorded in evidence/git-status-post.txt
- git diff --stat recorded in evidence/git-diff-stat-post.txt

## Validation
Script-only step. No runtime tests executed.
