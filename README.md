# dig-eco

## What this repo is
A minimal repository bootstrapped from zero using the Work Contract pattern.

## Local dev (Python 3.11)
Create a local venv and run tools from it:
- `python3.11 -m venv .venv`
- `.venv/bin/python -m pip install -U pip`
- `.venv/bin/python -m pip install ruff pytest`
- Tests: `.venv/bin/python -m pytest`
- Lint: `.venv/bin/python -m ruff check .`
- Format check: `.venv/bin/python -m ruff format --check .`

## Status & evidence
Audit artifacts live under `docs/status/audits/<STEP_ID>/`.
