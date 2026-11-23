#!/usr/bin/env python3
"""Turn a great_podflt run into IGS-named products: <AC>0<CMP><TYP>_<start>_<span>_05M_ORB.SP3 and _30S_CLK.CLK.

usage: great_product.py <run dir> --start YYYY-MM-DDTHH [--hours 24] [--ac GRT] [--cmp OPS] [--typ NRT] [--out DIR]
                        [--orb-step 300] [--clk-step 30]

Inputs (from the run directory): result/orbit.sp3 (or result/orbit_epochs.sp3, the per-epoch <flt> output, when the run
was stopped before the end) and clk_<yyyyddd>_epo (epoch-wise satellite clocks, RINEX clock 2 'AS' records).
The orbits are the forward-filter estimates; with a cold start the first ~12 h are not converged, so the product window
(--start, --hours) should cover the last part of the arc. SP3 clock values are filled from the clock file (microseconds;
999999.999999 where no clock).
"""
import argparse, datetime as dt, glob, os, sys

GPS0 = dt.datetime(1980, 1, 6)


def parse_sp3(path):
    epochs, cur = {}, None
    for ln in open(path):
        if ln.startswith('*'):
            f = ln[1:].split()
            sec = float(f[5])
            cur = dt.datetime(int(f[0]), int(f[1]), int(f[2]), int(f[3]), int(f[4])) + dt.timedelta(seconds=sec)
            epochs.setdefault(cur, {})
        elif ln.startswith('P') and cur is not None and len(ln) >= 46:
            epochs[cur][ln[1:4]] = (float(ln[4:18]), float(ln[18:32]), float(ln[32:46]))
    return epochs


def parse_clk(paths):
    clk = {}
    for p in paths:
        body = False
        for ln in open(p):
            if not body:
                body = 'END OF HEADER' in ln; continue
            if not ln.startswith('AS '): continue
            f = ln.split()
            sec = float(f[7])
            t = dt.datetime(int(f[2]), int(f[3]), int(f[4]), int(f[5]), int(f[6])) + dt.timedelta(seconds=sec)
            clk[(f[1], t)] = float(f[9])
    return clk


def name(A, start, kind):
    span = '%02dH' % A.hours if A.hours < 24 else '%02dD' % (A.hours // 24)
    smp = {'ORB': '%02dM' % (A.orb_step // 60), 'CLK': '%02dS' % A.clk_step if A.clk_step < 60 else '%02dM' % (A.clk_step // 60)}[kind]
    ext = 'SP3' if kind == 'ORB' else 'CLK'
    return '%s0%s%s_%s_%s_%s_%s.%s' % (A.ac, A.cmp, A.typ, start.strftime('%Y%j%H%M'), span, smp, kind, ext)


def write_sp3(path, A, orb, clk, epochs):
    sats = sorted({s for t in epochs for s in orb[t]}, key=lambda s: ('GRECJ'.index(s[0]) if s[0] in 'GRECJ' else 9, s))
    t0 = epochs[0]; d = t0 - GPS0
    wk, sow = d.days // 7, (d.days % 7) * 86400 + d.seconds
    mjd = (t0 - dt.datetime(1858, 11, 17)).days
    L = ['#dP%4d %2d %2d %2d %2d %11.8f %7d ORBIT IGS20 FIT %-4s' % (t0.year, t0.month, t0.day, t0.hour, t0.minute, t0.second, len(epochs), A.ac),
         '## %4d %15.8f %14.8f %5d %15.13f' % (wk, sow, A.orb_step, mjd, (t0.hour * 3600 + t0.minute * 60 + t0.second) / 86400.0)]
    nl = max(5, (len(sats) + 16) // 17)
    for l in range(nl):
        part = sats[l * 17:(l + 1) * 17]
        chunk = ''.join(part + ['  0'] * (17 - len(part)))
        L.append(('+  %3d   ' % len(sats) if l == 0 else '+        ') + chunk)
    for l in range(nl): L.append('++       ' + '  0' * 17)
    sysc = 'M' if len({s[0] for s in sats}) > 1 else {'G': 'G', 'R': 'R', 'E': 'E', 'C': 'C', 'J': 'J'}[sats[0][0]]
    L += ['%%c %s  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc' % sysc,
          '%c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc',
          '%f  1.2500000  1.025000000  0.00000000000  0.000000000000000',
          '%f  0.0000000  0.000000000  0.00000000000  0.000000000000000',
          '%i    0    0    0    0      0      0      0      0         0',
          '%i    0    0    0    0      0      0      0      0         0',
          '/* %s near-real-time orbits, GREAT-PODFLT filter (SRIF) on GEODNET stations' % A.ac,
          '/* forward-filter estimates; clocks from the same filter (reference: one station clock)']
    for t in epochs:
        L.append('*  %4d %2d %2d %2d %2d %11.8f' % (t.year, t.month, t.day, t.hour, t.minute, t.second + t.microsecond * 1e-6))
        for s in sats:
            if s not in orb[t]: continue
            x, y, z = orb[t][s]
            c = clk.get((s, t))
            L.append('P%s%14.6f%14.6f%14.6f%14.6f' % (s, x, y, z, c * 1e6 if c is not None else 999999.999999))
    L.append('EOF')
    with open(path, 'w') as f: f.write('\n'.join(L) + '\n')


def write_clk(path, A, clk, t0, t1):
    recs = sorted(((t, s, v) for (s, t), v in clk.items() if t0 <= t < t1 and
                   ((t - t0).total_seconds() % A.clk_step) == 0), key=lambda r: (r[0], r[1]))
    sats = sorted({r[1] for r in recs})
    L = ['     3.04           C                   M                   RINEX VERSION / TYPE',
         '%-20s%-20s%-20sPGM / RUN BY / DATE' % ('great_product.py', A.ac, dt.datetime.utcnow().strftime('%Y%m%d %H%M%S UTC')),
         '%-60sCOMMENT' % ('%s clocks from GREAT-PODFLT (GEODNET network)' % A.ac),
         '     1    AS                                                # / TYPES OF DATA',
         '%-3s  %-55sANALYSIS CENTER' % (A.ac, 'GEODNET GREAT-PODFLT'),
         '%6d%54sSAT / PRN COUNT' % (len(sats), ' ')]
    for i in range(0, len(sats), 15):
        L.append(''.join('%-4s' % s for s in sats[i:i + 15]).ljust(60) + 'PRN LIST')
    L.append('%60sEND OF HEADER' % ' ')
    for t, s, v in recs:
        L.append('AS %-4s %4d %02d %02d %02d %02d %9.6f  1   %19.12E' % (s, t.year, t.month, t.day, t.hour, t.minute,
                                                                     t.second + t.microsecond * 1e-6, v))
    with open(path, 'w') as f: f.write('\n'.join(L) + '\n')
    return len(recs)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run'); ap.add_argument('--start', required=True); ap.add_argument('--hours', type=int, default=24)
    ap.add_argument('--ac', default='GRT'); ap.add_argument('--cmp', default='OPS'); ap.add_argument('--typ', default='NRT')
    ap.add_argument('--out', default=None); ap.add_argument('--orb-step', type=int, default=300); ap.add_argument('--clk-step', type=int, default=30)
    A = ap.parse_args()
    start = dt.datetime.strptime(A.start, '%Y-%m-%dT%H'); end = start + dt.timedelta(hours=A.hours)
    res = os.path.join(A.run, 'result')
    sp3 = os.path.join(res, 'orbit.sp3')
    orb = parse_sp3(sp3) if os.path.exists(sp3) else {}
    if not orb:
        sp3 = os.path.join(res, 'orbit_epochs.sp3'); orb = parse_sp3(sp3)
    clk_paths = [os.path.join(res, 'satclk.clk')] if os.path.exists(os.path.join(res, 'satclk.clk')) else []
    clk_paths += sorted(glob.glob(os.path.join(A.run, 'clk_*_epo')))
    clk = parse_clk(clk_paths)
    epochs = [t for t in sorted(orb) if start <= t < end and (t - start).total_seconds() % A.orb_step == 0 and orb[t]]
    if not epochs: sys.exit('no orbit epochs in %s .. %s (%s)' % (start, end, sp3))
    out = A.out or os.path.join(A.run, 'products'); os.makedirs(out, exist_ok=True)
    p_orb = os.path.join(out, name(A, start, 'ORB')); write_sp3(p_orb, A, orb, clk, epochs)
    p_clk = os.path.join(out, name(A, start, 'CLK')); n = write_clk(p_clk, A, clk, start, end)
    print('%s: %d epochs (%s)\n%s: %d clock records' % (p_orb, len(epochs), os.path.basename(sp3), p_clk, n))


if __name__ == '__main__':
    main()
