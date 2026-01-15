# AGENTS

## Policy (Work Contract)
- Fail closed on ambiguity.
- Only modify files explicitly allowlisted for the current step.
- Capture evidence under docs/status/audits/<STEP_ID>/.
- Each step should be one bounded run -> evidence -> one commit.
- Allowlists govern tracked artifacts only; transient execution artifacts (e.g., `.venv/`) may be created during execution but must not be committed or added to allowlists.
