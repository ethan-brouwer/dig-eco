# dig-eco

## What this repo is
A minimal repository bootstrapped from zero using the Work Contract pattern.

## The workflow (Human -> ChatGPT -> Codex)
- Human defines intent, constraints, and acceptance.
- ChatGPT produces a canonical step packet (allowlist + commands + full file contents).
- Codex executes the packet, captures evidence, and produces one commit.

## Status & evidence
Audit artifacts live under `docs/status/audits/<STEP_ID>/`.

## How to continue
Use the Work Contract packet structure for every change:
intent -> allowlist -> execution -> evidence -> one commit.
