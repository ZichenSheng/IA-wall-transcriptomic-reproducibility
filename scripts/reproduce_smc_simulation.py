"""Optional Level B replay of seven Figure 6A M01 worlds with prepared input."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, nnls

from .common import ROOT, check, sha256


def run(path):
    expected = json.loads((ROOT / 'expected/smc.json').read_text())
    if sha256(path) != expected['sha256']:
        raise ValueError('Wrong prepared Figure 6A input. See DATA_SOURCES.md.')
    with np.load(path, allow_pickle=False) as a:
        reference, baseline, total = a['reference'], a['baseline'], a['total']
        if (baseline.shape != (300, 21) or reference.shape != (300, 6)
                or len(set(a['marker_genes'])) != 300 or not np.isfinite(reference).all()
                or list(a['classes']) != ['SMC', 'FB', 'EC', 'MAC', 'LYM', 'MONODC']):
            raise ValueError('Unexpected Figure 6 marker, specimen or reference frame')

        # Frozen E01 closure, pseudocount, weighted NNLS and coefficient normalization.
        weights = 1 / np.sqrt(reference.mean(axis=1) + 1)
        design = reference * weights[:, None]

        def shares(marker_counts, full_counts):
            if (full_counts <= 0).any():
                raise ValueError('Non-positive library size')
            marker_signal = 2 ** np.log2(marker_counts / full_counts * 1e6 + 0.5)
            normalized = marker_signal / marker_signal.sum(axis=0) * 1e6
            coefficients = np.array([nnls(design, y * weights)[0] for y in normalized.T])
            return coefficients[:, 0] / coefficients.sum(axis=1)

        original = shares(baseline, total)
        results = []
        for i, pair in enumerate(expected['pairs']):
            state = shares(baseline * a[f'state_{i}'][:, None], a[f'state_total_{i}'])
            delta = float((state - original).mean())
            pi, mass = a[f'pi_{i}'], a[f'smc_mass_{i}']

            def difference(factor):
                scaled = baseline * (1 + (factor - 1) * pi[:, None])
                scaled_total = total + (factor - 1) * mass
                return float((shares(scaled, scaled_total) - original).mean()) - delta

            if not difference(0.02) < 0 < difference(8):
                raise ValueError(f"Lost the frozen M01 bracket for pair {pair['pair_id']}")
            factor = brentq(difference, 0.02, 8, xtol=1e-13, rtol=1e-13, maxiter=200)
            check(factor, pair['matched_factor'], 5e-8)
            check(delta, pair['delta_share'], 5e-9)
            check(difference(factor), 0, 1e-12)
            results.append({'pair_id': pair['pair_id'], 'category': pair['category'],
                            'equivalent_scaling': factor, 'mean_smc_share_change': delta})
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prepared_input', type=Path)
    args = parser.parse_args()
    print(json.dumps({'status': 'PASS', 'level': 'B', 'worlds': run(args.prepared_input)}, indent=2))
