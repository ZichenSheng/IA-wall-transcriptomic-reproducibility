"""Figure 3A observed statistic from prepared Welch vectors and Level C inference checks."""
import numpy as np

from ._d1_kernel import LES, RUPT, dstat, fast_rin
from .common import check, read


def run():
    contrasts = [f'R:{c}' for c in RUPT] + [f'L:{c}' for c in LES]
    frame = read('rupture_welch.tsv.gz', float_precision='round_trip').set_index('gene', verify_integrity=True)
    if list(frame.columns) != contrasts or not frame.index.is_monotonic_increasing:
        raise ValueError('Unexpected gene or contrast order')
    if not np.isfinite(frame.to_numpy()).all():
        raise ValueError('Non-finite Welch vector')
    check(len(frame), 15648)
    counts = read('rupture_counts.tsv').set_index('cohort', verify_integrity=True)
    groups = counts.loc[RUPT, ['RIA', 'UIA']].to_numpy()
    if groups.tolist() != [[21,21],[11,8],[8,6],[8,5]]:
        raise ValueError('Unexpected rupture specimen group counts')

    fit = dstat({key: fast_rin(frame[key].to_numpy()) for key in contrasts})
    without = float(np.mean([excess for pair, excess in fit['D'].items() if 'GSE54083' not in pair]))
    frozen = read('rupture_frozen.tsv').set_index('variant')
    draws = read('rupture_draws.tsv.gz')
    result = {'n_genes': len(frame), 'n_rupture_comparison_specimens': int(groups.sum())}
    for label, key, obs in [('primary','PRIMARY',fit['Dbar']),
                            ('drop','S4_DROP_RUPTURE_GSE54083',without)]:
        row = frozen.loc[key]
        null = draws[f'{label}_null'].dropna().to_numpy()
        boot = draws[f'{label}_boot'].dropna().to_numpy()
        check(len(null), row.n_null)
        check(len(boot), 2000 if label == 'primary' else 1000)
        check(obs, row.Dbar, 5e-10)
        tail = (1 + int((null >= obs).sum())) / (len(null) + 1)
        low, high = np.quantile(boot, [0.025, 0.975])
        check(tail, row.p, 5e-10)
        check(low, row.lo, 5e-10)
        check(high, row.hi, 5e-10)
        result[label] = {'observed_level_B': obs, 'archived_draws_level_C':
                         {'p': tail, 'ci': [float(low), float(high)], 'n_null': len(null)}}
    return result
