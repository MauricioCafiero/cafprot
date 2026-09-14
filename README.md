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

## Python 3.14 and PROPKA (handled automatically)

**PROPKA 3.5.1 — the newest release that exists — crashes on Python 3.14.** It
reads instance `__annotations__` at runtime to type its parameter file, and
PEP 649 made annotations lazy:

```
AttributeError: 'Parameters' object has no attribute '__annotations__'.
Did you mean: '__annotate_func__'?
```

The trap is that pdb2pqr still *succeeds* without PROPKA — it accepts
`--with-ph`, exits 0, and does no pKa prediction at all, so a 3.14 run can look
perfectly fine while ignoring the pH entirely.

`cafprot` fixes this rather than documenting around it. The annotations still
exist and `annotationlib` can materialise them; reattaching them to the *class*
doesn't help, because 3.14 redirects that assignment into
`__annotations_cache__` where instance lookup won't find it — but putting the
dict into each instance's `__dict__` does work. That shim is applied inside the
pdb2pqr subprocess, where PROPKA actually runs.

**Result: pH-aware receptor protonation works on both interpreters**, and on
3.14 the output is byte-identical to 3.11's native PROPKA at pH 1, 7 and 13.

| | 3.11.15 | 3.14.5 |
|---|---|---|
| `normalize_smiles` | works | works |
| `protonate_receptor`, hydrogens added | works | works |
| `protonate_receptor`, **pH applied** | yes (native) | yes (via shim) |

If PROPKA ever becomes genuinely unrunnable, `propka_available()` returns False
and `protonate_receptor()` warns loudly that the pH was not applied, instead of
returning a file that quietly pretends otherwise.

## Tests

```sh
python tests/test_cafprot.py     # or: pytest tests/
```

11 tests, no pytest required, no network, and no dependency on any other repo —
the only input is the 4-residue fragment in `tests/mini_peptide.pdb`. Verified
11/11 on both 3.11.15 and 3.14.5.

`test_receptor_ph_is_actually_applied` asserts that pH 1 and pH 13 give
different protonation, so an interpreter where PROPKA silently does nothing
fails the suite rather than passing it.

One test pins a behaviour rather than endorsing it: phenol (pKa ≈ 10) comes back
as `[O-]c1ccccc1` at pH 7.4, because Dimorphite-DL's phenol rule is aggressive.
That is the library's call, not this module's; it is pinned so a version bump
that changes it gets noticed.
