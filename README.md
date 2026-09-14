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
cafprot.protonation_states("Oc1ccccc1")                # -> ['Oc1ccccc1', '[O-]c1ccccc1']
cafprot.protonate_receptor("receptor.pdb")             # -> 'receptor_protonated.pdb'
cafprot.propka_available()                             # -> True / False
```

| | Default | Notes |
|---|---|---|
| `normalize_smiles(smiles, ph, canonicalize_tautomer=False)` | `ph=PH_LIGAND` = **7.4** | Dimorphite-DL; blood plasma; one state |
| `protonation_states(smiles, ph, precision=1.0)` | `ph=PH_LIGAND` = **7.4** | every plausible state; see the disclaimer |
| `protonate_receptor(pdb_path, ph, overwrite=False)` | `ph=PH_RECEPTOR` = **7.0** | pdb2pqr + PROPKA; pdb2pqr's own default |

The two pH defaults live as constants at the top of `cafprot.py`, so the
family-wide convention is one edit rather than a decision re-made in every CLI.

All three **fall back to returning their input unchanged, with a warning,
rather than raising** — one bad molecule should never abort a batch. Receptor
output is cached as `<stem>_protonated.pdb` beside the input and reused unless
`overwrite=True`.

## pH range, and how far to trust it

**Any pH from 0 to 14 works** — the physiological defaults are only defaults.
Values outside that range are rejected by pdb2pqr, and `protonate_receptor()`
then warns and returns the input unchanged rather than crashing.

Measured with `normalize_smiles`:

| pH | aspirin (COOH, pKa ≈ 3.5) | ethylamine (pKa ≈ 10.7) | imidazole |
|---|---|---|---|
| 0–2 | neutral | `[NH3+]` | `[nH+]` |
| 4–6 | `[O-]` | `[NH3+]` | `[nH+]` → neutral |
| 7.4 | `[O-]` | `[NH3+]` | `[n-]` |
| 9–14 | `[O-]` | neutral | `[n-]` |

and with `protonate_receptor` on the ASP/GLU test fragment: 25 hydrogens at
pH 0–2, 24 from pH 4 up — the carboxylates titrating right at pKa ≈ 4. It is
flat above 4 only because that fragment has nothing else ionisable.

**Disclaimer: `normalize_smiles` returns ONE state, and buying that single
answer costs you the uncertainty Dimorphite-DL was built to express.**

Dimorphite-DL associates each of its 38 ionizable moieties with a pKa *range*,
`[µ − nσ, µ + nσ]`, not a point value — deliberately, because it treats each
ionizable site independently and the ranges are how it absorbs that
approximation. When the range overlaps the requested pH it emits *both* the
protonated and deprotonated forms, meaning "either is plausible."

To return a single string, `normalize_smiles` asks for `precision=0` (n = 0),
which collapses each range to its bare mean. For chemically heterogeneous
classes that mean is a poor estimate for any specific molecule:

| Molecule | Class | µ used | true pKa | result at 7.4 |
|---|---|---|---|---|
| phenol | `Phenol` = `[c,n,o:1]-[O:2]-[H]` | 7.07 (σ 3.28) | ≈ 10 | ionised |
| imidazole | `[n:1]-[H]` | 7.17 (σ 2.95) | ≈ 14.5 | `[n-]` |
| ethylamine | `[C:1]-[NX3+0:2]` | 8.16 (σ 2.52) | 10.7 | neutral by pH 9 |

Those SMARTS are broad — `Phenol` also matches N–OH and O–OH and electron-poor
phenols near pKa 7; the aromatic N–H class spans tetrazole (≈ 4.9) to pyrrole
(≈ 17) — so σ of 2.5–3.3 log units is the library flagging a heterogeneous
class, and `precision=0` discards that flag. **This is `cafprot`'s choice, not
a defect in Dimorphite-DL.**

Well-defined groups are fine: aspirin flips between pH 2 and 4, exactly where
pKa 3.5 predicts. For a ligand carrying a phenol, an azole N–H or another
group in a broad class, call `protonation_states()` to see every form the
library considers plausible before trusting the single answer:

```python
cafprot.protonation_states("Oc1ccccc1", ph=7.4)   # ['Oc1ccccc1', '[O-]c1ccccc1']
```

The receptor side is unaffected — it goes through PROPKA, which titrated the
test carboxylates correctly at pKa ≈ 4.

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

12 tests, no pytest required, no network, and no dependency on any other repo —
the only input is the 4-residue fragment in `tests/mini_peptide.pdb`. Verified
12/12 on both 3.11.15 and 3.14.5.

`test_receptor_ph_is_actually_applied` asserts that pH 1 and pH 13 give
different protonation, so an interpreter where PROPKA silently does nothing
fails the suite rather than passing it.

One test pins a behaviour rather than endorsing it: phenol (pKa ≈ 10) comes back
as `[O-]c1ccccc1` at pH 7.4, because the single-answer API collapses the Phenol
class to its 7.07 mean. Pinned so a version bump that shifts it gets noticed,
and paired with a test asserting `protonation_states()` still surfaces both
forms.

## License

[MIT](LICENSE)
