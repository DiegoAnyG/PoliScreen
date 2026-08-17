# Contributing

Issues, questions and pull requests are welcome.

## Reporting something

Open an issue with:

- what you ran — the command, or the steps in the interface,
- what you expected and what happened,
- the output of `poliscreen info`, and the **Versions** block of an exported Methods file if the
  problem is about results.

For a wrong or surprising *result*, attach the project's `run.json` and the ranking. Results depend
on the receptor, the control and the residue roles, so the configuration matters more than the
error message.

## Pull requests

1. `python -m pytest tests/ -q` passes. CI runs it too.
2. A change that only shows up at runtime brings a test with it — Streamlit widget keys, packaging,
   Windows paths and residue numbering all have one already, because each of them broke once.
3. Code, comments, docstrings and file names in English. The Spanish that is there is deliberate:
   the interface catalogue, the legacy on-disk names, and the user-facing CSV aliases.
4. Comments say *why*, not what. The repository is a record of decisions and the reasons they were
   made; a fix without its reason gets undone by the next person who finds it odd.

## What will not be merged

Anything that bends a result toward a conclusion. Parameters may be changed for determinism or
convergence, and a scoring change is welcome **with the measurement that justifies it** — on
enough cases to be believable. During development a change was reverted here after one favourable
case turned into 46 lost interactions across thirty.

## Support

Maintained by one person, alongside a doctorate. Issues are read; a fix arrives when it arrives.
Security or data-loss reports are answered first.
