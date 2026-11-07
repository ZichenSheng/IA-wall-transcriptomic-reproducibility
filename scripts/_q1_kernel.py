"""Frozen Q1 arithmetic from p0_q1_core and p0_dependency."""
import numpy as np
from scipy import stats

K = 200
B = 2000
BBOOT = 10000


def pool_mu(y, v):   # R2.1 DL pooling, mu and se only (verbatim arithmetic)
    y, v = np.asarray(y, float), np.asarray(v, float); k = len(y)
    if k == 1: return y[0], np.sqrt(v[0]), 0.0
    w = 1 / v; mu_f = (w * y).sum() / w.sum(); Q = (w * (y - mu_f) ** 2).sum()
    c = w.sum() - (w ** 2).sum() / w.sum(); tau2 = max(0.0, (Q - (k - 1)) / c)
    ws = 1 / (v + tau2); mu = (ws * y).sum() / ws.sum(); I2 = max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0
    return mu, np.sqrt(1 / ws.sum()), I2


def pool_many(Y, V):
    """vectorised DL over rows (genes x cohorts); identical arithmetic to pool_mu."""
    k = Y.shape[1]
    if k == 1: return Y[:, 0]
    w = 1 / V; mu_f = (w * Y).sum(1) / w.sum(1); Q = (w * (Y - mu_f[:, None]) ** 2).sum(1)
    c = w.sum(1) - (w ** 2).sum(1) / w.sum(1); tau2 = np.maximum(0.0, (Q - (k - 1)) / c)
    ws = 1 / (V + tau2[:, None]); return (ws * Y).sum(1) / ws.sum(1)


def pctile(obs, y): return float(((y < obs).sum() + 0.5 * (y == obs).sum()) / len(y))


def null_pcts(obs, y):
    """matched percentile of each pool member j against {obs} U y\\{y_j} (K comparators)."""
    s = np.sort(y); out = np.empty(len(y))
    for j, v in enumerate(y):
        less = np.searchsorted(s, v, 'left') - 0; eq = np.searchsorted(s, v, 'right') - np.searchsorted(s, v, 'left') - 1   # exclude itself
        less += (obs < v); eq += (obs == v); out[j] = (less + 0.5 * eq) / len(y)
    return out


def global_null(NP, idx, seed):
    """verbatim c02 global_test null draw sequence."""
    g = np.random.default_rng(seed); return np.array([np.mean([NP[i][g.integers(K)] for i in idx]) for _ in range(B)])


def cluster_stats(d, cl, seed, B=BBOOT):
    cl = np.asarray(cl); labs, inv = np.unique(cl, return_inverse=True); G = len(labs); N = len(d)
    Dg = np.bincount(inv, weights=d, minlength=G); ng = np.bincount(inv, minlength=G).astype(float); delta = d.mean()
    se = np.sqrt(G / (G - 1) * ((Dg - delta * ng) ** 2).sum()) / N; se_iid = d.std(ddof=1) / np.sqrt(N)
    rng = np.random.default_rng(seed); bs = np.empty(B)
    for b in range(B):
        s = rng.integers(G, size=G); bs[b] = Dg[s].sum() / ng[s].sum()
    lo, hi = np.quantile(bs, [0.025, 0.975]); t = delta / se
    return dict(n_units=N, n_clusters=G, largest_cluster=int(ng.max()), delta=round(delta, 4), theta=round(delta + 0.5, 4), ci95_lo=round(lo, 4), ci95_hi=round(hi, 4),
                se_cluster=round(se, 5), se_iid=round(se_iid, 5), design_effect=round((se / se_iid) ** 2, 2), t_cluster=round(t, 2), df=G - 1,
                p_one_sided_t=float(stats.t.sf(t, G - 1)), boot_frac_le_0=float((bs <= 0).mean()), seed=seed, B=B, passes=bool(lo > 0))
