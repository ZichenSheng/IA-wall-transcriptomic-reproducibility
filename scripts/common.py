"""Shared file identity and result checks."""
import hashlib
import json
import math
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'


def sha256(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def validate_inputs(root=ROOT):
    manifest = json.loads((root / 'expected/input_sha256.json').read_text())
    observed = {str(p.relative_to(root)) for p in (root / 'data').rglob('*') if p.is_file()}
    if observed != set(manifest):
        raise ValueError('Unexpected or missing input files')
    for rel, digest in manifest.items():
        path = root / rel
        if path.is_symlink() or sha256(path) != digest:
            raise ValueError(f'Input identity failed for {rel}')


def read(name, **kwargs):
    return pd.read_csv(DATA / name, sep='\t', **kwargs)


def check(value, expected, tolerance=0):
    if not math.isclose(value, expected, rel_tol=0, abs_tol=tolerance):
        raise ValueError(f'Observed {value}, expected {expected} (absolute tolerance {tolerance})')
