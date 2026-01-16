# Evidence: 20260115T235052Z-mrds-gee-csv-fix

## Scope
Fix MRDS GEE script for CSV geometry/commodity parsing and add El Salvador boundary.

## Pre-state
- git status recorded in evidence/git-status-pre.txt

## Changes
- Updated `mrds_gee_map.js` to build geometry from X/Y and parse commod1 from description.
- Added El Salvador boundary layer.
- Kept commodity legend with counts.

## Post-state
- git status recorded in evidence/git-status-post.txt
- git diff --stat recorded in evidence/git-diff-stat-post.txt

## Validation
Script-only step. No runtime tests executed.
