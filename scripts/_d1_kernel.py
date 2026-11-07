"""Frozen D1 arithmetic from d1lib and d1_engine."""
from itertools import combinations
import numpy as np
from scipy.special import ndtri

RUPT = ['GSE122897', 'GSE13353', 'GSE15629', 'GSE54083']
LES = ['GSE122897', 'GSE54083', 'GSE15629', 'GSE26969']
PAIRS = list(combinations(RUPT, 2))
SUB3 = ['GSE122897', 'GSE13353', 'GSE54083']
PAIRS3 = list(combinations(SUB3, 2))


def tvec(X, a, b):
    """Gene-wise Welch t of mean(a)-mean(b); X genes x samples ndarray; a,b column index lists."""
    A = X[:, a]; B = X[:, b]
    se = np.sqrt(A.var(1, ddof=1) / A.shape[1] + B.var(1, ddof=1) / B.shape[1])
    with np.errstate(divide='ignore', invalid='ignore'):
        return (A.mean(1) - B.mean(1)) / se


def zs(v):
    v = np.asarray(v, float)
    return (v - v.mean()) / v.std()


def one_factor(Rm, n_iter=500, tol=1e-12):
    """Minres one-factor fit on a correlation matrix (off-diagonal least squares). Returns loadings."""
    K = Rm.shape[0]
    # start from first principal axis
    w, V = np.linalg.eigh(Rm)
    l = V[:, -1] * np.sqrt(max(w[-1], 1e-6))
    if l.sum() < 0:
        l = -l
    off = ~np.eye(K, dtype=bool)
    for _ in range(n_iter):
        l_old = l.copy()
        for k in range(K):
            o = [j for j in range(K) if j != k]
            num = sum(Rm[k, j] * l[j] for j in o); den = sum(l[j] ** 2 for j in o)
            l[k] = num / den if den > 0 else 0.0
        if np.max(np.abs(l - l_old)) < tol:
            break
    return np.clip(l, -0.999, 0.999)


def fast_rin(v):
    n = len(v)
    r = np.empty(n); r[np.argsort(v, kind='mergesort')] = np.arange(1, n + 1)
    z = ndtri((r - 0.375) / (n + 0.25))
    return (z - z.mean()) / z.std()


def cmat(Z, names):
    M = np.array([Z[n] for n in names])
    return (M @ M.T) / M.shape[1]


def dstat(Z, lesion=None, pair_specific=True, loadings=None, alpha_mode=False):
    """Z: dict of standardized vectors with keys R:<c>, L:<k>. Returns dict of statistics.
    lesion: list of lesion dataset ids to use (default LES present in Z)."""
    les = [k for k in (lesion or LES) if f'L:{k}' in Z]
    Ln = [f'L:{k}' for k in les]
    RL = {}
    n = len(next(iter(Z.values())))
    for x in RUPT:
        for k in les:
            RL[(x, k)] = float(Z[f'R:{x}'] @ Z[f'L:{k}'] / n)
    LL = cmat(Z, Ln)
    if loadings is None:
        ell = dict(zip(les, one_factor(LL))) if len(les) >= 3 else None
    else:
        ell = loadings
    res = dict(ell=ell, LL=LL, les=les)
    D, pred, rr, mm = {}, {}, {}, {}
    for i, j in PAIRS:
        rr[(i, j)] = float(Z[f'R:{i}'] @ Z[f'R:{j}'] / n)
        if pair_specific:
            Ls = [k for k in les if k not in (i, j)]
        else:
            Ls = None
        def m_of(x, Lset):
            Lset = [k for k in (Lset if Lset is not None else les) if k != x]
            if not Lset:
                return np.nan
            if alpha_mode or ell is None:
                # composite with Cronbach-alpha disattenuation
                comp = zs(np.mean([Z[f'L:{k}'] for k in Lset], axis=0))
                K = len(Lset)
                if K == 1:
                    return np.nan
                sub = cmat(Z, [f'L:{k}' for k in Lset]); rbar = sub[~np.eye(K, dtype=bool)].mean()
                alpha = K * rbar / (1 + (K - 1) * rbar)
                if alpha <= 0:
                    return np.nan
                return float(Z[f'R:{x}'] @ comp / n) / np.sqrt(alpha)
            num = sum(ell[k] * RL[(x, k)] for k in Lset); den = sum(ell[k] ** 2 for k in Lset)
            return num / den
        mi = m_of(i, Ls); mj = m_of(j, Ls)
        mm[(i, j)] = (mi, mj)
        pred[(i, j)] = mi * mj
        D[(i, j)] = rr[(i, j)] - pred[(i, j)]
    res.update(D=D, pred=pred, rr=rr, mpair=mm, RL=RL)
    vals = np.array([D[p] for p in PAIRS])
    res['Dbar'] = float(np.nanmean(vals))
    res['Dbar3'] = float(np.nanmean([D[p] for p in PAIRS3]))
    res['Pi'] = float(np.nansum([pred[p] for p in PAIRS]) / np.nansum([rr[p] for p in PAIRS]))
    res['Pi3'] = float(np.nansum([pred[p] for p in PAIRS3]) / np.nansum([rr[p] for p in PAIRS3]))
    res['rbar'] = float(np.mean([rr[p] for p in PAIRS])); res['rbar3'] = float(np.mean([rr[p] for p in PAIRS3]))
    # cohort-level loadings on lambda using all lesions disjoint from the cohort
    mc = {}
    for x in RUPT:
        Ls = [k for k in les if k != x]
        if ell is not None and not alpha_mode:
            mc[x] = sum(ell[k] * RL[(x, k)] for k in Ls) / sum(ell[k] ** 2 for k in Ls)
        else:
            mc[x] = np.nan
    res['m_c'] = mc
    # extra-factor loadings within the 3-cohort subset
    a, b, c = SUB3
    Dab, Dac, Dbc = D[(a, b)], D[(a, c)], D[(b, c)]
    p2 = {}
    for x, (u, v, w) in {a: (Dab, Dac, Dbc), b: (Dab, Dbc, Dac), c: (Dac, Dbc, Dab)}.items():
        p2[x] = u * v / w if w != 0 else np.nan
    res['p2'] = p2
    return res
