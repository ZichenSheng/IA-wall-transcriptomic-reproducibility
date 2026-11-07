"""Figure 2A and primary co-expression cluster bootstrap in Figure 2B."""
import numpy as np
from scipy.stats import false_discovery_control

from ._q1_kernel import B, K, cluster_stats, global_null, null_pcts, pctile, pool_many, pool_mu
from .common import DATA, check, read


def run():
    claims = read('q1/claims.tsv', keep_default_na=False)
    frozen = read('q1/frozen_claim_results.tsv')
    check(len(claims), 239)
    check(claims.gene.nunique(), 218)
    endpoint_counts = claims.endpoint.value_counts().to_dict()
    if endpoint_counts != {'IA_VS_CONTROL': 124, 'RIA_VS_UIA': 64, 'UIA_VS_CONTROL': 33, 'RIA_VS_CONTROL': 18}:
        raise ValueError('Unexpected endpoint counts')
    if claims.unit_id.duplicated().any() or not claims.unit_id.equals(frozen.unit_id):
        raise ValueError('Claim ID order or uniqueness failed')
    if not (claims.unit_id == claims.gene + '|' + claims.endpoint).all():
        raise ValueError('Claim gene and endpoint mismatch')

    backgrounds = {}
    for path in sorted((DATA / 'q1').glob('BG_*.tsv.gz')):
        cohort, endpoint = path.name[3:-7].split('_', 1)
        table = read(f'q1/{path.name}').set_index('HGNC_ID', verify_integrity=True)
        if not np.isfinite(table.to_numpy()).all() or not (table.v > 0).all():
            raise ValueError(f'Invalid effect background for {cohort}/{endpoint}')
        backgrounds[cohort, endpoint] = table

    percentiles, null_percentiles, signed_g, p = [], [], [], []
    for row in claims.itertuples(index=False):
        cohorts = row.cohorts_unexposed_used.split(';')
        measured = set(row.cohorts_measured.split(';'))
        exposed = set(row.cohorts_exposed.split(';'))
        if cohorts != sorted(measured - exposed) or row.direction not in {'UP', 'DOWN'}:
            raise ValueError(f'Invalid claim direction or cohort assignment for {row.unit_id}')
        tables = [backgrounds[cohort, row.endpoint] for cohort in cohorts]
        gene = row.hgnc_id
        common = tables[0].index
        for table in tables[1:]:
            common = common.intersection(table.index)
        if gene not in common or len(common) <= K:
            raise ValueError(f'Incomplete matching universe for {row.unit_id}')
        # Frozen expr_pool rule. The first sorted unexposed cohort defines distance.
        covariates = tables[0].loc[common, ['pct_mean', 'pct_sd']].to_numpy(float)
        target = tables[0].loc[gene, ['pct_mean', 'pct_sd']].to_numpy(float)
        distance = np.sqrt(((covariates - target) ** 2).sum(axis=1))
        distance[common.get_loc(gene)] = np.inf
        matched = common.to_numpy()[np.argpartition(distance, K)[:K]]
        sign = 1 if row.direction == 'UP' else -1
        effect, _, _ = pool_mu([sign * table.at[gene, 'g'] for table in tables],
                               [table.at[gene, 'v'] for table in tables])
        y = np.column_stack([table.loc[matched, 'g'].to_numpy(float) for table in tables])
        v = np.column_stack([table.loc[matched, 'v'].to_numpy(float) for table in tables])
        controls = sign * pool_many(y, v)
        percentiles.append(pctile(effect, controls))
        null_percentiles.append(null_pcts(effect, controls))
        signed_g.append(effect)
        p.append((1 + int((controls >= effect).sum())) / (K + 1))

    percentiles = np.array(percentiles)
    np.testing.assert_array_equal(np.round(percentiles, 4), frozen.expr_matched_percentile)
    np.testing.assert_allclose(np.round(signed_g, 4), frozen.signed_pooled_g, atol=1e-9, rtol=0)
    q = false_discovery_control(p)
    np.testing.assert_allclose(q, frozen.expr_matched_p_BH_descriptive, atol=5.1e-7, rtol=0)
    obs = float(np.round(percentiles, 4).mean())
    null = global_null(null_percentiles, np.arange(len(claims)), 20261101)
    draws_ge_obs = int((null >= obs).sum())
    delta = percentiles - np.array([n.mean() for n in null_percentiles])
    boot = cluster_stats(delta, claims.cluster.to_numpy(), 2026102701)
    check(round(obs, 4), 0.7378)
    check(round(float(null.mean()), 4), 0.4991)
    check(draws_ge_obs, 0)
    check(boot['n_clusters'], 31)
    check(boot['ci95_lo'], 0.1245)
    check(boot['ci95_hi'], 0.2876)
    check(int((q < 0.05).sum()), 0)
    sign_reversals = int((np.array(signed_g) < 0).sum())
    check(sign_reversals, 43)
    return {'n_claims': len(claims), 'n_genes': claims.gene.nunique(), 'mean_percentile': obs,
            'matched_null_mean': float(null.mean()), 'matched_null_p': (1 + draws_ge_obs) / (B + 1),
            'minimum_BH_q': float(q.min()), 'sign_reversals': sign_reversals, 'endpoint_counts': endpoint_counts,
            'cluster_bootstrap': {'clusters': boot['n_clusters'], 'reps': boot['B'],
                                  'theta': boot['theta'], 'ci': [round(boot[k] + 0.5, 4)
                                                                for k in ('ci95_lo','ci95_hi')]}}
