#!/usr/bin/env python3
"""Ocean tide loading BLQ file for GREAT from FES2004 spherical-harmonic synthesis (FES2004S1.dat up to degree/order 100,
load Love numbers Love_load_cm.dat, CM frame).

usage: make_blq.py <sites.snx> <models dir with FES2004S1.dat + Love_load_cm.dat> <out.blq> [YYYY-MM-DD]

The synthesis is linear in each tidal constituent, so the BLQ amplitude/phase of the 11 BLQ waves (M2 S2 N2 K2 K1 O1 P1 Q1 Mf Mm
Ssa) follow exactly from evaluating the constituent's loading-potential coefficients at argument 0 and 90 deg. Phases are given
for GREAT's argument convention (t_gtideIERS::load_ocean: Schwiderski arguments, displacement = A cos(chi - phi), rows radial,
west, south); the constant offset between the Doodson astronomical argument and GREAT's Schwiderski argument of each wave is
measured numerically.
GREAT applies no nodal corrections, so the nodal factor f and angle u of the date (default: today; f/u follow the 18.6-year
lunar node, so regenerate the file every few months) are folded into the amplitudes and phases.
Not representable in a BLQ (and so dropped): the extra constituents of the synthesis (Sa, S1, 2N2, M4, Mtm, Msq, Om1, Om2).

Records carry '$$ ... lon/lat: <lon> <lat> <h>' because GREAT matches BLQ records by position (10 km), not by name.
The output holds station positions: write it outside the public repository.
"""
import math, sys
import numpy as np

KMAX = 100
GM, AE, KF = 3.986004415e14, 6378137.0, 1.0 / 298.25641153
D2R = math.pi / 180
WAVES = [('M2', 255555), ('S2', 273555), ('N2', 245655), ('K2', 275555), ('K1', 165555), ('O1', 145555), ('P1', 163555),
         ('Q1', 135655), ('MF', 75555), ('MM', 65455), ('SSA', 57555)]
BIAS = {55565: 180, 125755: -90, 127555: -90, 135655: -90, 137455: -90, 145545: -90, 145555: -90, 155555: -90, 155655: 90,
        157455: 90, 162556: -90, 163555: -90, 163655: 90, 164556: 90, 165555: 90, 166554: 90, 167555: 90, 175455: 90, 185555: 90,
        253755: 180, 255545: 180, 265455: 180, 263655: 180, 275554: 180, 355555: 180}


def load_tables(d):
    fes = []
    with open(d + '/FES2004S1.dat') as f:
        for i, ln in enumerate(f):
            if i < 3: continue
            p = ln.split()
            if len(p) < 12: continue
            n, m = int(p[2]), int(p[3])
            if n > KMAX or n < 1: continue
            fes.append((int(round(float(p[0]) * 1e3)), n, m, float(p[8]), float(p[9]), float(p[10]), float(p[11])))
    love = [(0.0, 0.0, 0.0)]
    with open(d + '/Love_load_cm.dat') as f:
        for i, ln in enumerate(f):
            if i < 6: continue
            p = ln.replace('D', 'E').replace('d', 'E').split()
            if len(p) < 4 or len(love) > KMAX: break
            love.append((float(p[1]), float(p[2]), float(p[3])))
    return fes, love


def nodal(mjd, dood):
    """nodal factor f and phase u (deg) of a constituent (FES2004/IERS nodal tide formulation)"""
    T = (mjd + 2400000.5 - 2451545.0) / 36525.0
    N = (((2.22222e-6 * T + 2.0708e-3) * T - 1934.136261) * T + 125.04452) % 360.0
    sn, cn, s2n, c2n = math.sin(N * D2R), math.cos(N * D2R), math.sin(2 * N * D2R), math.cos(2 * N * D2R)
    F, U = [1.0] * 23, [0.0] * 23
    F[1] = F[2] = 1.009 + 0.187 * cn - 0.015 * c2n
    F[3] = F[4] = 1.006 + 0.115 * cn - 0.009 * c2n
    F[5] = F[6] = 1.000 - 0.037 * cn
    F[7] = F[8] = 1.024 + 0.286 * cn + 0.008 * c2n
    F[9] = F[10] = F[11] = math.hypot(1 + 0.189 * cn - 0.0058 * c2n, 0.189 * sn - 0.0058 * s2n)
    F[12] = math.hypot(1 + 0.185 * cn, 0.185 * sn); F[13] = math.hypot(1 + 0.198 * cn, 0.198 * sn)
    F[14] = math.hypot(1 + 0.640 * cn + 0.134 * c2n, 0.640 * sn + 0.134 * s2n)
    F[15] = F[16] = F[17] = F[18] = math.hypot(1 - 0.0373 * cn, 0.0373 * sn)
    F[19] = F[3]; F[20] = F[6] * math.sqrt(F[6]); F[21] = F[6] ** 2; F[22] = F[6] ** 3
    U[1] = U[2] = 10.8 * sn - 1.3 * s2n
    U[3] = U[4] = -8.9 * sn + 0.7 * s2n
    U[5] = U[6] = -2.1 * sn
    U[7] = U[8] = -17.7 * sn + 0.7 * s2n
    U[9] = U[10] = U[11] = math.atan2(0.189 * sn - 0.0058 * s2n, 1 + 0.189 * cn - 0.0058 * s2n) / D2R
    U[12] = math.atan2(0.185 * sn, 1 + 0.185 * cn) / D2R; U[13] = math.atan2(-0.198 * sn, 1 + 0.198 * cn) / D2R
    U[14] = math.atan2(-0.640 * sn - 0.134 * s2n, 1 + 0.640 * cn + 0.134 * c2n) / D2R
    U[15] = U[16] = U[17] = U[18] = math.atan2(-0.0373 * sn, 1 - 0.0373 * cn) / D2R
    U[19] = U[3]; U[20] = 1.5 * U[6]; U[21] = 2 * U[6]; U[22] = 3 * U[6]
    dd = [0, 135655, 145555, 163555, 165555, 245655, 255555, 273555, 275555, 125755, 127555, 137455, 155655, 175455, 185555,
          235755, 237555, 247455, 265455, 164556, 355555, 455555, 655555]
    for i in range(1, 23):
        if dd[i] == dood: return F[i], U[i]
    return 1.0, 0.0


def wave_coefs(fes, dood):
    """loading-potential coefficients (cnm, snm) of one constituent at argument 0 and at 90 deg"""
    fk = 4 * math.pi * 6.67428e-11 * 1.025e3 / 9.7803278 * 1e-2
    size = (KMAX + 2) * (KMAX + 2) + 2
    out = []
    for th in (0.0, math.pi / 2):
        c, s = np.zeros(size), np.zeros(size)
        for (d, n, m, cpa, cpp, cma, cmp) in fes:
            if d != dood: continue
            k = n * (n + 1) // 2 + m
            cp, sp = cpa * math.sin(cpp * D2R) / (2 * n + 1), cpa * math.cos(cpp * D2R) / (2 * n + 1)
            cn, sn = cma * math.sin(cmp * D2R) / (2 * n + 1), cma * math.cos(cmp * D2R) / (2 * n + 1)
            c[k] += (cp + cn) * math.cos(th) + (sp + sn) * math.sin(th)
            s[k] -= (cp - cn) * math.sin(th) - (sp - sn) * math.cos(th)
        out.append((c * fk, s * fk))
    return out


def legendre(t):
    N = (KMAX + 2) * (KMAX + 2) + 2
    pnm, dp = np.zeros(N), np.zeros(N)
    u = math.sqrt(1 - t * t)
    pnm[1], pnm[2], pnm[3] = 1.0, math.sqrt(3) * t, math.sqrt(3) * u
    for n in range(1, KMAX + 1):
        a = math.sqrt(2 * n + 1) / math.sqrt(2 * n - 1)
        b = math.sqrt(2 * (n - 1) * (2 * n + 1)) / math.sqrt((2 * n - 1) * n)
        kk, k1, k2 = n * (n + 1) // 2 + 1, n * (n - 1) // 2 + 1, n * (n - 1) // 2 + 2
        pnm[kk] = a * t * pnm[k1] - b * u / 2 * pnm[k2]
        for m in range(1, n + 1):
            kk, k1, k2, k3 = n * (n + 1) // 2 + m + 1, n * (n - 1) // 2 + m + 1, n * (n - 1) // 2 + m + 2, n * (n - 1) // 2 + m
            c = a / n * math.sqrt(n * n - m * m)
            d = a / n / 2 * math.sqrt((n - m) * (n - m - 1))
            e = a / n / 2 * math.sqrt((n + m) * (n + m - 1)) * (math.sqrt(2) if m == 1 else 1)
            pnm[kk] = c * t * pnm[k1] - d * u * pnm[k2] + e * u * pnm[k3]
    dp[2], dp[3] = -math.sqrt(3) * u, math.sqrt(3) * t
    for n in range(2, KMAX + 1):
        kk, k1 = n * (n + 1) // 2 + 1, n * (n + 1) // 2 + 2
        dp[kk] = -math.sqrt(n * (n + 1) / 2) * pnm[k1]
        kk, k1, k2 = n * (n + 1) // 2 + 2, n * (n + 1) // 2 + 1, n * (n + 1) // 2 + 3
        dp[kk] = math.sqrt(2 * n) * math.sqrt(n + 1) / 2 * pnm[k1] - math.sqrt(n - 1) * math.sqrt(n + 2) / 2 * pnm[k2]
        for m in range(2, n + 1):
            kk, k1, k2 = n * (n + 1) // 2 + m + 1, n * (n + 1) // 2 + m, n * (n + 1) // 2 + m + 2
            dp[kk] = math.sqrt(n + m) * math.sqrt(n - m + 1) / 2 * pnm[k1] - math.sqrt(n - m) * math.sqrt(n + m + 1) / 2 * pnm[k2]
    return pnm, dp


def geodetic(x, y, z):   # WGS84
    a, f = 6378137.0, 1 / 298.257223563
    e2 = f * (2 - f); lon = math.atan2(y, x); p = math.hypot(x, y); lat = math.atan2(z, p * (1 - e2)); h = 0.0
    for _ in range(10):
        N = a / math.sqrt(1 - e2 * math.sin(lat) ** 2); h = p / math.cos(lat) - N
        lat = math.atan2(z, p * (1 - e2 * N / (N + h)))
    return lat, lon, h


class Site:
    def __init__(self, xyz, love):
        lat, lon, h = geodetic(*xyz)
        e2 = (2 - KF) * KF; sl = math.sin(lat); N = AE / math.sqrt(1 - e2 * sl * sl)
        x, y = (N + h) * math.cos(lat) * math.cos(lon), (N + h) * math.cos(lat) * math.sin(lon)
        z = (N * (1 - e2) + h) * sl
        self.lat, self.lon, self.h = lat, lon, h
        self.rr = math.sqrt(x * x + y * y + z * z)
        self.rlat = math.atan(z / math.hypot(x, y)); self.rlon = lon
        self.gr = 9.7803253359 * (1 + 0.00193185265241 * sl * sl) / math.sqrt(1 - 0.00669437999013 * sl * sl) - 3.086e-6 * h
        pnm, dp = legendre(math.sin(self.rlat))
        # per packed index j: geometry factors for east, north, up
        ns, ms, js = [], [], []
        for n in range(1, KMAX + 1):
            for m in range(n + 1):
                ns.append(n); ms.append(m); js.append(n * (n + 1) // 2 + m)
        ns, ms, js = np.array(ns), np.array(ms), np.array(js)
        sc = (AE / self.rr) ** ns
        dh = np.array([love[n][0] for n in ns]); dl = np.array([love[n][1] for n in ns])
        cl, sl_ = np.cos(self.rlon * ms), np.sin(self.rlon * ms)
        self.js = js
        self.fe = (ms * pnm[js + 1] * dl * sc, cl, sl_)
        self.fn = dp[js + 1] * dl * sc
        self.fu = pnm[js + 1] * dh * sc
        self.k = GM / self.rr / self.gr

    def enu(self, cnm, snm):
        c, s = cnm[self.js], snm[self.js]
        me, cl, sl = self.fe
        te = np.sum(me * (c * sl - s * cl)); tn = np.sum((c * cl + s * sl) * self.fn); tu = np.sum((c * cl + s * sl) * self.fu)
        return np.array([-te * self.k / math.cos(self.rlat), -tn * self.k, tu * self.k])


# ---- arguments: Doodson astronomical (theta) and GREAT Schwiderski (chi), radians, t in MJD (UTC)
def theta(dood, mjd):
    t = (mjd - 51544.5) / 36525.0; tt = (mjd - 51544.0) - math.floor(mjd - 51544.0)
    a = [0] * 7
    a[2] = 218.31664563 + 481267.88119575 * t - 0.001466388889 * t * t - 0.000000074112 * t ** 3 - 0.000000153389 * t ** 4
    a[3] = 280.4664501606 + 36000.769748805556 * t + 0.000303222222 * t * t - 0.000001905501 * t ** 3 - 0.000000065361 * t ** 4
    a[4] = 83.35324312 + 4069.01363525 * t - 0.010321722222 * t * t - 0.000014417168 * t ** 3 + 0.0000000526333 * t ** 4
    a[5] = 234.95544499 + 1934.136261972222 * t - 0.002075611111 * t * t - 0.000000213944 * t ** 3 + 0.000000164972 * t ** 4
    a[6] = 282.9373409806 + 1.719457666668 * t + 0.000456888889 * t * t - 0.000001943279 * t ** 3 - 0.0000000033444 * t ** 4
    a[1] = 360.0 * tt - a[2] + a[3]
    td, bs, th = dood, 100000, 0.0
    for j in range(1, 7):
        d = td // bs; td -= bs * d; bs //= 10
        if j > 1: d -= 5
        th += d * a[j]
    return (th + BIAS.get(dood, 0.0)) * D2R


SPEED = [1.40519e-4, 1.45444e-4, 1.37880e-4, 1.45842e-4, 0.72921e-4, 0.67598e-4, 0.72523e-4, 0.64959e-4, 0.053234e-4, 0.026392e-4,
         0.003982e-4]
ANGFAC = [[2, 0, 2, 2, 1, 1, -1, 1, 0, 0, 2], [-2, 0, -3, 0, 0, -2, 0, -3, 2, 1, 0], [0, 0, 1, 0, 0, 0, 0, 1, 0, -1, 0],
          [0, 0, 0, 0, 0.25, -0.25, -0.25, -0.25, 0, 0, 0]]


def chi(i, mjd):
    import datetime as dt
    t = dt.datetime(1858, 11, 17) + dt.timedelta(days=mjd)
    doy = t.timetuple().tm_yday
    capd = doy + 365 * (t.year - 1975) + (t.year - 1973) // 4
    capt = (27392.500528 + 1.000000035 * capd) / 36525.0
    h0 = (279.69668 + (36000.768930485 + 3.03e-4 * capt) * capt) * D2R
    s0 = (((1.9e-6 * capt - 0.001133) * capt + 481267.88314137) * capt + 270.434358) * D2R
    p0 = (((-1.2e-5 * capt - 0.010325) * capt + 4069.0340329577) * capt + 334.329653) * D2R
    sod = t.hour * 3600 + t.minute * 60 + t.second + t.microsecond * 1e-6
    return SPEED[i] * sod + ANGFAC[0][i] * h0 + ANGFAC[1][i] * s0 + ANGFAC[2][i] * p0 + ANGFAC[3][i] * 2 * math.pi


def wrap(a):
    return (a + math.pi) % (2 * math.pi) - math.pi


def offsets():
    """chi - theta per wave (constant; checked over two years)"""
    out = []
    for i, (nm, d) in enumerate(WAVES):
        v = [wrap(chi(i, m) - theta(d, m)) for m in np.arange(61000.0, 61730.0, 7.37)]
        ref = v[0]; dev = max(abs(wrap(x - ref)) for x in v)
        if dev > 2 * D2R: sys.exit('argument offset of %s not constant (%.2f deg)' % (nm, dev / D2R))
        out.append(math.atan2(np.mean(np.sin(np.array(v))), np.mean(np.cos(np.array(v)))))
    return out


def read_snx(path):
    xyz, on = {}, False
    for ln in open(path):
        if ln.startswith('+SOLUTION/ESTIMATE'): on = True; continue
        if ln.startswith('-SOLUTION/ESTIMATE'): break
        if on and ln[7:10] == 'STA':
            f = ln.split(); xyz.setdefault(f[2], [0, 0, 0])['XYZ'.index(f[1][3])] = float(f[8])
    return xyz


def main():
    snx, mdir, out = sys.argv[1:4]
    import datetime as dt
    day = dt.datetime.strptime(sys.argv[4], '%Y-%m-%d') if len(sys.argv) > 4 else dt.datetime.utcnow()
    mjd = (day - dt.datetime(1858, 11, 17)).total_seconds() / 86400.0
    fu = [nodal(mjd, d) for _, d in WAVES]
    fes, love = load_tables(mdir)
    coefs = [wave_coefs(fes, d) for _, d in WAVES]
    off = offsets()
    sites = read_snx(snx)
    L = ['$$ Ocean loading displacement, FES2004 harmonic synthesis (CM frame), make_blq.py',
         '$$ phases for GREAT t_gtideIERS::load_ocean (Schwiderski arguments); nodal f/u of %s included' % day.strftime('%Y-%m-%d'),
         '$$ COLUMN ORDER:  M2  S2  N2  K2  K1  O1  P1  Q1  MF  MM SSA',
         '$$ ROW ORDER: AMPLITUDES (m) RADIAL, TANGENTL EW (west), TANGENTL NS (south); PHASES (deg) same order',
         '$$ END HEADER']
    for name in sorted(sites):
        s = Site(sites[name], love)
        amp, ph = np.zeros((3, 11)), np.zeros((3, 11))
        for i, ((c0, s0), (c9, s9)) in enumerate(coefs):
            dc, ds = s.enu(c0, s0), s.enu(c9, s9)          # d(t) = dc cos(theta) + ds sin(theta), ENU
            for r, (a, b) in enumerate(((dc[2], ds[2]), (-dc[0], -ds[0]), (-dc[1], -ds[1]))):   # up, west, south
                R = math.hypot(a, b); alpha = math.atan2(b, a)  # = R cos(theta - alpha) = R cos(chi - (alpha + off))
                f, u = fu[i]   # f R cos(theta + u - alpha) = f R cos(chi - (alpha - u + off))
                amp[r, i] = f * R; ph[r, i] = wrap(alpha - u * D2R + off[i]) / D2R
        lon = s.lon / D2R
        L.append('$$')
        L.append('  %s' % name)
        L.append('$$ %s FES2004 synthesis lon/lat: %10.4f %10.4f %10.3f' % (name, lon, s.lat / D2R, s.h))
        for r in range(3): L.append(' ' + ''.join(' %.5f' % v for v in amp[r]))
        for r in range(3): L.append(' ' + ''.join(' %6.1f' % v for v in ph[r]))
    L.append('$$ END TABLE')
    with open(out, 'w') as f: f.write('\n'.join(L) + '\n')
    print('%s: %d sites' % (out, len(sites)))


if __name__ == '__main__':
    main()
