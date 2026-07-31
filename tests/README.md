# Tests

```bash
conda activate poliscreen && pytest tests/ -q
```

Takes a few seconds. Needs no Vina, ADCP, gnina or network: it covers the logic,
not the engines.

## What it covers

- **`test_core.py`** — peptide chemistry (head-to-tail cyclization, terminus
  protection, net charge), Windows path translation, resource allocation by
  memory and flexibility, peptide recognition from structure, and in-memory
  export.
- **`test_interfaz.py`** — renders every step and every ligand mode looking for
  exceptions, with an empty folder and with data, and checks that state survives
  switching steps.

Each test documents in its docstring **the real failure that motivated it**. A
test without that context ends up deleted when it gets in the way.

## What it does not cover

Anything depending on external binaries and a specific target: docking itself,
PLIP, AGFR preparation and ADMET prediction. That needs a real run; these tests
only guarantee that nothing around them is broken.

## When adding one

Write the test that reproduces the failure first, see it red, then fix. Three
quarters of this application's failures have been silent — a plausible but wrong
result, not an exception — so check the **value**, not just that it does not
crash.
