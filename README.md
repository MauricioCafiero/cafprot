# cafprot

pH-aware protonation for ligands (SMILES) and receptors (PDB), in one place.

Everything in the surrounding project family currently adds hydrogens by naive
rules — RDKit's `AddHs()`, which honours whatever charges the input SMILES
happened to carry, or `obabel -h`, which is valence-based rather than pKa-based.
Neither asks what charge state a molecule or a residue actually has at
physiological pH. This answers that.

**Not a package.** There is nothing to `pip install`, no `pyproject.toml`, no
version to pin. Copy `cafprot.py` next to the code that needs it, or put this
repo on `sys.path`.

**Not integrated with anything.** By design — no other repo imports this yet.

## API

```python
import cafprot

cafprot.normalize_smiles("CC(=O)Oc1ccccc1C(=O)O")     # -> 'CC(=O)Oc1ccccc1C(=O)[O-]'
cafprot.normalize_smiles("CCN")                        # -> 'CC[NH3+]'
cafprot.protonate_receptor("receptor.pdb")             # -> 'receptor_protonated.pdb'
cafprot.propka_available()                             # -> True / False
```

| | Default | Notes |
|---|---|---|
| `normalize_smiles(smiles, ph, canonicalize_tautomer=False)` | `ph=PH_LIGAND` = **7.4** | Dimorphite-DL; blood plasma |
| `protonate_receptor(pdb_path, ph, overwrite=False)` | `ph=PH_RECEPTOR` = **7.0** | pdb2pqr + PROPKA; pdb2pqr's own default |

The two pH defaults live as constants at the top of `cafprot.py`, so the
family-wide convention is one edit rather than a decision re-made in every CLI.

Both functions **fall back to returning their input unchanged, with a warning,
rather than raising** — one bad molecule should never abort a batch. Receptor
output is cached as `<stem>_protonated.pdb` beside the input and reused unless
`overwrite=True`.

## Setup

```sh
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Python 3.14 caveat (measured, not assumed)

**PROPKA 3.5.1 crashes on Python 3.14.** It reads instance `__annotations__` at
runtime, which PEP 649 changed:

```
AttributeError: 'Parameters' object has no attribute '__annotations__'.
Did you mean: '__annotate_func__'?
```

The trap is that pdb2pqr still *succeeds* on 3.14 — it accepts `--with-ph` and
exits 0 while doing no pKa prediction at all, so a run can look fine while
ignoring the pH entirely. `protonate_receptor()` detects this via
`propka_available()` and warns explicitly that the pH was not applied, instead
of returning a file that quietly pretends otherwise.

**Use Python 3.11 for pH-aware receptor work.** The ligand side is unaffected
and works on both.

| | 3.11.15 | 3.14.5 |
|---|---|---|
| `normalize_smiles` | works | works |
| `protonate_receptor`, hydrogens added | works | works |
| `protonate_receptor`, pH applied | **yes** | **no** (warns) |

## Tests

```sh
python tests/test_cafprot.py     # or: pytest tests/
```

10 tests, no pytest required, no network, and no dependency on any other repo —
the only input is the 4-residue fragment in `tests/mini_peptide.pdb`. Verified
10/10 on both 3.11.15 and 3.14.5.

One test pins a behaviour rather than endorsing it: phenol (pKa ≈ 10) comes back
as `[O-]c1ccccc1` at pH 7.4, because Dimorphite-DL's phenol rule is aggressive.
That is the library's call, not this module's; it is pinned so a version bump
that changes it gets noticed.
