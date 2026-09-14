"""pH-aware protonation for ligands (SMILES) and receptors (PDB).

Everything in the surrounding project family currently adds hydrogens by
naive rules -- RDKit's AddHs() (which honours whatever charges the input
SMILES happened to carry) or `obabel -h` (valence-based, not pKa-based).
Neither asks what charge state a molecule or a residue actually has at
physiological pH. This module answers that question in one place:

    normalize_smiles(smiles, ph)   ligand side, via Dimorphite-DL
    protonate_receptor(pdb, ph)    receptor side, via pdb2pqr + PROPKA

Not a package -- there is nothing to pip install. Copy cafprot.py next to
the code that needs it, or put this repo on sys.path.

Both functions fall back to returning their input unchanged, with a warning,
rather than raising: one bad molecule should never abort a batch.

PYTHON 3.14 CAVEAT (measured, not assumed): PROPKA 3.5.1 crashes on 3.14 --
it reads instance __annotations__ at runtime, which PEP 649 changed. pdb2pqr
itself still runs, and will even accept --with-ph while silently doing no
pKa prediction at all, so a 3.14 run can look successful while ignoring the
pH entirely. protonate_receptor() detects this and says so loudly. For real
pH-aware receptor prep, use Python 3.11.
"""
import os
import shutil
import subprocess
import sys
import warnings

PH_LIGAND = 7.4      # blood plasma
PH_RECEPTOR = 7.0    # pdb2pqr's own default, closer to intracellular


def normalize_smiles(smiles, ph=PH_LIGAND, canonicalize_tautomer=False):
    """The dominant protonation state of `smiles` at `ph`, as canonical SMILES.

    Returns `smiles` unchanged (with a warning) if Dimorphite-DL is missing or
    cannot handle the input. Set canonicalize_tautomer=True to additionally
    collapse tautomers via RDKit -- a related but distinct problem, off by
    default because it changes structures protonation alone would leave alone.
    """
    from rdkit import Chem

    try:
        from dimorphite_dl import protonate_smiles
        # precision=0 asks for the single dominant state; the default of 1.0
        # returns every variant within a +/-1 pKa window instead.
        states = protonate_smiles(smiles, ph_min=ph, ph_max=ph, precision=0.0)
    except Exception as exc:
        warnings.warn(f"protonation failed for {smiles!r} ({exc}); using input unchanged")
        return smiles

    if not states:
        warnings.warn(f"no protonation state returned for {smiles!r}; using input unchanged")
        return smiles

    mol = Chem.MolFromSmiles(states[0])
    if mol is None:
        warnings.warn(f"unparseable protonated form {states[0]!r} for {smiles!r}; using input unchanged")
        return smiles

    if canonicalize_tautomer:
        from rdkit.Chem.MolStandardize import rdMolStandardize
        mol = rdMolStandardize.TautomerEnumerator().Canonicalize(mol)

    return Chem.MolToSmiles(mol)


def propka_available():
    """Whether PROPKA can actually run here -- i.e. whether a requested pH will
    mean anything. False on Python 3.14; see the module docstring.
    """
    try:
        import propka.parameters
        propka.parameters.Parameters().__annotations__
        return True
    except Exception:
        return False


def protonate_receptor(pdb_path, ph=PH_RECEPTOR, overwrite=False):
    """A properly protonated copy of `pdb_path`, written as
    `<stem>_protonated.pdb` next to the input, whose path is returned.

    Cached: an existing output is reused unless overwrite=True. Returns the
    input path unchanged (with a warning) if pdb2pqr is missing or fails.

    Without PROPKA the pH CANNOT be applied; rather than return a file that
    silently ignores it, this warns explicitly and still returns the
    hydrogen-added structure, which remains better than obabel's valence-only
    guess.
    """
    out_pdb = f"{os.path.splitext(pdb_path)[0]}_protonated.pdb"
    if os.path.exists(out_pdb) and not overwrite:
        return out_pdb

    # Console scripts live beside the interpreter, which matters when a venv's
    # python is invoked directly (venv/bin/python ...) without activation --
    # then the venv's bin/ is not on PATH and a plain which() finds nothing.
    exe = None
    for path in (None, os.path.dirname(sys.executable)):
        exe = shutil.which("pdb2pqr", path=path) or shutil.which("pdb2pqr30", path=path)
        if exe:
            break
    if exe is None:
        warnings.warn("pdb2pqr not found on PATH; using the unprotonated input")
        return pdb_path

    cmd = [exe, "--ff=AMBER", f"--pdb-output={out_pdb}"]
    if propka_available():
        cmd += [f"--with-ph={ph}", "--titration-state-method=propka"]
    else:
        warnings.warn(
            f"PROPKA unavailable under this interpreter, so pH {ph} was NOT applied to "
            f"{pdb_path}; hydrogens were added without pKa prediction. Use Python 3.11 "
            f"for pH-aware receptor protonation."
        )
    cmd += [pdb_path, f"{os.path.splitext(pdb_path)[0]}_protonated.pqr"]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 or not os.path.exists(out_pdb):
        warnings.warn(
            f"pdb2pqr failed on {pdb_path} ({result.stderr.strip().splitlines()[-1:] or 'no stderr'}); "
            f"using the unprotonated input"
        )
        return pdb_path
    return out_pdb
