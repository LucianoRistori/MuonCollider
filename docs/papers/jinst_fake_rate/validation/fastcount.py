"""
fastcount.py - exact enumeration of ALL hit combinations passing a chi2 cut
(straight lines, two views), without trying all prod(n_i) combinations.

Planes are added one at a time. For a partial combination, the chi2 of the
fit to its hits can only grow when a plane is added, so a partial
combination already above the cut can never pass: dropping it loses
nothing. Adding a hit at (u, v) on the next plane raises chi2 by exactly
    [(u - u_hat)^2 + (v - v_hat)^2] / (sigma^2 + V)
(u_hat, v_hat: the partial fit's prediction there, V: its variance), so
only hits inside a disc of radius sqrt((cut - chi2_partial)(sigma^2 + V))
around the prediction are tried. The result is identical to fitting every
combination (checked against brute force in validate_fastcount.py).
"""
import numpy as np


class Plane:
    def __init__(self, Y, W, sigma):
        self.Y, self.W, self.sigma = float(Y), float(W), float(sigma)


def search(planes, hits, cut):
    """planes: list of Plane (the set S, in order); hits: list of (u, v)
    arrays in mm, one pair per plane. Returns an int array (n_pass, |S|) of
    hit indices of every combination with chi2 < cut, and their chi2."""
    m = len(planes)
    # sort each plane's hits by u for window searches
    order = [np.argsort(h[0]) for h in hits]
    us = [hits[i][0][order[i]] for i in range(m)]
    vs = [hits[i][1][order[i]] for i in range(m)]
    # first two planes: every pair (chi2 = 0 for two points)
    n0, n1 = len(us[0]), len(us[1])
    i0 = np.repeat(np.arange(n0), n1)
    i1 = np.tile(np.arange(n1), n0)
    idx = [i0, i1]
    w0, w1 = 1 / planes[0].sigma ** 2, 1 / planes[1].sigma ** 2
    Y0, Y1 = planes[0].Y, planes[1].Y
    # running weighted sums, per candidate, for both views
    S0 = np.full(len(i0), w0 + w1)
    S1 = np.full(len(i0), w0 * Y0 + w1 * Y1)
    S2 = np.full(len(i0), w0 * Y0 ** 2 + w1 * Y1 ** 2)
    Su = w0 * us[0][i0] + w1 * us[1][i1]
    SuY = w0 * Y0 * us[0][i0] + w1 * Y1 * us[1][i1]
    Sv = w0 * vs[0][i0] + w1 * vs[1][i1]
    SvY = w0 * Y0 * vs[0][i0] + w1 * Y1 * vs[1][i1]
    chi2 = np.zeros(len(i0))
    for k in range(2, m):
        Yk, sk = planes[k].Y, planes[k].sigma
        wk = 1 / sk ** 2
        D = S0 * S2 - S1 ** 2
        bu = (S0 * SuY - S1 * Su) / D
        au = (Su - bu * S1) / S0
        bv = (S0 * SvY - S1 * Sv) / D
        av = (Sv - bv * S1) / S0
        uh, vh = au + bu * Yk, av + bv * Yk
        V = (S2 - 2 * Yk * S1 + Yk ** 2 * S0) / D
        denom = sk ** 2 + V
        rad = np.sqrt(np.maximum(cut - chi2, 0) * denom)
        lo = np.searchsorted(us[k], uh - rad, "left")
        hi = np.searchsorted(us[k], uh + rad, "right")
        cnt = hi - lo
        if cnt.sum() == 0:
            return np.zeros((0, m), int), np.zeros(0)
        cand = np.repeat(np.arange(len(lo)), cnt)
        offs = np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt)
        j = lo[cand] + offs
        du = us[k][j] - uh[cand]
        dv = vs[k][j] - vh[cand]
        new = chi2[cand] + (du ** 2 + dv ** 2) / denom[cand]
        keep = new < cut
        cand, j, new = cand[keep], j[keep], new[keep]
        idx = [a[cand] for a in idx] + [j]
        S0 = S0[cand] + wk
        S1 = S1[cand] + wk * Yk
        S2 = S2[cand] + wk * Yk ** 2
        Su = Su[cand] + wk * us[k][j]
        SuY = SuY[cand] + wk * Yk * us[k][j]
        Sv = Sv[cand] + wk * vs[k][j]
        SvY = SvY[cand] + wk * Yk * vs[k][j]
        chi2 = new
    # back to the caller's hit numbering
    out = np.stack([order[i][idx[i]] for i in range(m)], 1) if len(chi2) else np.zeros((0, m), int)
    return out, chi2
