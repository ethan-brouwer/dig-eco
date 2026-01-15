# Evidence: 20260115T222113Z-python-scaffold

## Scope
Apply preset python311-venv-ruff-pytest: create minimal Python scaffold, validate with pytest and ruff, and record audit evidence.

## Tooling
This step uses a local virtualenv (.venv) for ruff/pytest (not committed).
- Python: evidence/python311-version.txt
- Installed packages: evidence/pip-freeze.txt

## Pre-state
- git status: evidence/git-status-pre.txt

## Validation
- `.venv/bin/python -m pytest -q`: evidence/pytest.txt
- `.venv/bin/python -m ruff check .`: evidence/ruff-check.txt
- `.venv/bin/python -m ruff format --check .`: evidence/ruff-format-check.txt

## Post-state
- git status: evidence/git-status-post.txt
- git diff --stat: evidence/git-diff-stat-post.txt
