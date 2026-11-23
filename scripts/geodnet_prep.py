#!/usr/bin/env python3
"""Prepare a GREAT-PODFLT run on GEODNET stations: RTCM3 logs -> RINEX, station SINEX, ERP, alias map, podflt.xml.

usage: geodnet_prep.py --stations <csv> --beg YYYY-MM-DDTHH --end YYYY-MM-DDTHH --work <dir> [options]

Everything this script writes is PRIVATE (station IDs, coordinates) and goes below --work, which must be outside the
public repository. The repository holds only this script and the XML template (config/podflt_geodnet.xml.in).

Inputs
  --stations   CSV with header 'name,x,y,z[,...]' (12-character GEODNET IDs, ECEF metres). These coordinates are held fixed
               (written to the SINEX with sigma 1 mm).
  --rtcm-root  hourly 30 s RTCM3 logs <root>/<YYYY-MM-DD>-30/<ID>_<YYYY>_<MM>_<DD>_30_<hh>_<hh+1>.bin (hourly binary
               logs; broadcast ephemerides from BCEP00BKG0 logs in the same directories)
  --erp        IGS ERP version 2 file covering the arc (e.g. from IERS Bulletin A); daily records are taken
  --sys-dir    directory with the GREAT system files (jpleph_de405, EIGEN6S4, ocean_tide_FES2004, ...);
               --atx / --satinfo override the ANTEX / satinfo files (use the ANTEX release of the reference products)

Station names: GREAT reads 4-character site codes only (SINEX, BLQ and RINEX parsers). Each GEODNET ID gets a stable
4-character alias (<work>/alias_map.csv, created on first use and only appended to, never re-ordered). The RINEX marker
name and the file names use the alias.

Outputs (below --work)
  alias_map.csv, sites.snx, erp_great, leap_seconds, obs/<YYYY>/<DDD>/<alias><DDD>0.<YY>o, nav/<YYYY>/<DDD>/brdm<DDD>0.<YY>p,
  run_<tag>/podflt.xml (+ result/ and logger/)
"""
import argparse, concurrent.futures as cf, csv, datetime as dt, os, shutil, string, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
ALPHA = string.digits + string.ascii_uppercase


def parse_t(s):
    return dt.datetime.strptime(s, '%Y-%m-%dT%H') if 'T' in s else dt.datetime.strptime(s, '%Y-%m-%d')


def alias_for(i):
    # 'G' + 3 base-36 digits: 46656 aliases; the leading G keeps them distinct from IGS-style 4-letter names
    return 'G' + ALPHA[(i // 1296) % 36] + ALPHA[(i // 36) % 36] + ALPHA[i % 36]


def load_alias(path, names):
    amap = {}
    if os.path.exists(path):
        for r in csv.DictReader(open(path)):
            amap[r['id']] = r['alias']
    used = set(amap.values())
    i, new = 0, []
    for n in names:
        if n in amap: continue
        while alias_for(i) in used: i += 1
        amap[n] = alias_for(i); used.add(amap[n]); new.append(n)
    if new or not os.path.exists(path):
        with open(path, 'w') as f:
            f.write('id,alias\n')
            for k, v in amap.items(): f.write('%s,%s\n' % (k, v))
        os.chmod(path, 0o600)
    return amap


def hours_of_day(day, beg, end):
    d0 = dt.datetime(day.year, day.month, day.day)
    return [h for h in range(24) if d0 + dt.timedelta(hours=h + 1) > beg and d0 + dt.timedelta(hours=h) < end]


def concat_logs(root, name, day, hours, out):
    n = 0
    with open(out, 'wb') as o:
        for h in hours:
            p = os.path.join(root, day.strftime('%Y-%m-%d') + '-30',
                             '%s_%s_30_%03d_%03d.bin' % (name, day.strftime('%Y_%m_%d'), h, h + 1))
            if os.path.exists(p) and os.path.getsize(p) > 0:
                with open(p, 'rb') as i: shutil.copyfileobj(i, o); n += 1
    return n


def strip_rinex_antenna(path):
    # GEODNET coordinates refer to the GPS L1+L2 iono-free APC: antenna type ADVNULLANTENNA (not in the ANTEX, so no
    # PCO/PCV is applied) and zero H/E/N delta, whatever antenna the RTCM stream announces.
    with open(path, 'rb') as f:
        data = f.read()
    lines = data.splitlines(keepends=True)
    for i, ln in enumerate(lines):
        if len(ln) >= 60:
            lbl = ln[60:80]
            if b'ANT # / TYPE' in lbl:
                eol = b'\r\n' if ln.endswith(b'\r\n') else b'\n'
                lines[i] = b'%-20s%-20s%20sANT # / TYPE        ' % (b'0', b'ADVNULLANTENNA', b'') + eol
            elif b'ANTENNA: DELTA H/E/N' in lbl:
                eol = b'\r\n' if ln.endswith(b'\r\n') else b'\n'
                lines[i] = b'        0.0000        0.0000        0.0000                  ANTENNA: DELTA H/E/N' + eol
        if b'END OF HEADER' in ln:
            break
    tmp = path + '.hdrtmp'
    with open(tmp, 'wb') as f:
        f.writelines(lines)
    os.replace(tmp, path)


def to_rinex_obs(A, name, alias, day, beg, end):
    hours = hours_of_day(day, beg, end)
    odir = os.path.join(A.work, 'obs', day.strftime('%Y'), day.strftime('%j'))
    obs = os.path.join(odir, '%s%s0.%so' % (alias.lower(), day.strftime('%j'), day.strftime('%y')))
    if os.path.exists(obs) and not A.force: return 'skip'
    tmp = os.path.join(A.work, 'tmp', '%s_%s.bin' % (alias, day.strftime('%j')))
    if concat_logs(A.rtcm_root, name, day, hours, tmp) == 0:
        os.remove(tmp)
        return 'kept' if os.path.exists(obs) else 'nodata'   # logs already cleaned up: keep the earlier conversion
    d = day.strftime('%Y/%m/%d')
    subprocess.run([A.convbin, '-r', 'rtcm3', '-tr', d, '12:00:00', '-ts', d, '%02d:00:00' % hours[0], '-te', d,
                    '%02d:59:30' % hours[-1], '-ti', '30', '-v', '3.04', '-hm', alias, '-ha', '0/ADVNULLANTENNA', '-od', '-os', '-o', obs, tmp],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    os.remove(tmp)
    if not os.path.exists(obs): return 'fail'
    strip_rinex_antenna(obs)
    n = 0
    with open(obs, 'rb') as f:
        for ln in f:
            if ln.startswith(b'> '): n += 1
    if n < A.min_epochs:
        os.remove(obs); return 'short'
    return 'ok'


def to_rinex_nav(A, day, beg, end):
    ndir = os.path.join(A.work, 'nav', day.strftime('%Y'), day.strftime('%j')); os.makedirs(ndir, exist_ok=True)
    nav = os.path.join(ndir, 'brdm%s0.%sp' % (day.strftime('%j'), day.strftime('%y')))
    if os.path.exists(nav) and not A.force: return nav
    # whole day (the nav file is selected by day in GREAT): BCEP00BKG0 plus the extra relay streams given
    tmp = os.path.join(A.work, 'tmp', 'nav_%s.bin' % day.strftime('%j'))
    with open(tmp, 'wb') as o:
        for src in ['BCEP00BKG0'] + [s for s in A.nav_streams.split(',') if s]:
            part = tmp + '.' + src
            if concat_logs(A.rtcm_root, src, day, range(24), part):
                with open(part, 'rb') as i: shutil.copyfileobj(i, o)
            os.remove(part)
    d = day.strftime('%Y/%m/%d')
    out_obs = os.path.join(A.work, 'tmp', 'nav_%s.obs' % day.strftime('%j'))
    new = nav + '.new'
    if os.path.getsize(tmp) > 0:
        subprocess.run([A.convbin, '-r', 'rtcm3', '-tr', d, '12:00:00', '-v', '3.04', '-o', out_obs, '-n', new, tmp],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for p in (tmp, out_obs):
        if os.path.exists(p): os.remove(p)
    # replace only with a non-empty conversion: the RTCM logs are cleaned up after a few days, keep the earlier file then
    if os.path.exists(new) and os.path.getsize(new) > 10000: os.replace(new, nav)
    elif os.path.exists(new): os.remove(new)
    if not os.path.exists(nav): print('warning: no broadcast navigation for %s' % day.strftime('%Y-%m-%d'))
    return nav


def write_sinex(path, sites, beg, end):
    def snx_t(t):
        return '%02d:%03d:%05d' % (t.year % 100, t.timetuple().tm_yday, t.hour * 3600 + t.minute * 60 + t.second)
    b, e = snx_t(beg - dt.timedelta(days=2)), snx_t(end + dt.timedelta(days=2))
    now = snx_t(dt.datetime.utcnow())
    L = ['%%=SNX 2.02 GRT %s GRT %s %s P %05d 0 S' % (now, b, e, 3 * len(sites)),
         '+FILE/REFERENCE', ' DESCRIPTION        GEODNET station coordinates for GREAT-PODFLT (geodnet_prep.py)', '-FILE/REFERENCE',
         '+SITE/ID', '*CODE PT __DOMES__ T _STATION DESCRIPTION__ APPROX_LON_ APPROX_LAT_ _APP_H_']
    for a, _ in sites: L.append(' %-4s  A %9s P %-22s %3d %02d %04.1f %3d %02d %04.1f %7.1f' % (a, '', 'GEODNET', 0, 0, 0, 0, 0, 0, 0))
    L += ['-SITE/ID', '+SOLUTION/EPOCHS', '*CODE PT SOLN T _DATA_START_ __DATA_END__ _MEAN_EPOCH_']
    for a, _ in sites: L.append(' %-4s  A    1 P %s %s %s' % (a, b, e, b))
    L += ['-SOLUTION/EPOCHS', '+SOLUTION/ESTIMATE', '*INDEX _TYPE_ CODE PT SOLN _REF_EPOCH__ UNIT S ___ESTIMATED_VALUE___ __STD_DEV__']
    k = 0
    for a, xyz in sites:
        for comp, v in zip(('STAX', 'STAY', 'STAZ'), xyz):
            k += 1
            L.append(' %5d %-6s %-4s  A    1 %s m    2 %21.14E %11.5E' % (k, comp, a, b, v, 0.001))
    L += ['-SOLUTION/ESTIMATE', '%ENDSNX']
    with open(path, 'w') as f: f.write('\n'.join(L) + '\n')
    os.chmod(path, 0o600)


def write_erp(src, dst, beg, end, leap=37):
    # IGS ERP v2 (MJD, X, Y [1e-6"], UT1-UTC [1e-7 s]) -> GREAT daily table, full UT1 (mode 'UT1', not 'UT1R': the
    # tidal terms are removed / restored by gtrs2crs)
    recs = {}
    for ln in open(src):
        f = ln.split()
        try: mjd = float(f[0])
        except (ValueError, IndexError): continue
        if abs(mjd - round(mjd)) < 1e-6 and len(f) >= 4:
            try:
                recs[int(round(mjd))] = (float(f[1]) * 1e-6, float(f[2]) * 1e-6, float(f[3]) * 1e-7 - leap)
            except ValueError:
                continue
    m0 = int((beg - dt.datetime(1858, 11, 17)).days) - 2
    m1 = int((end - dt.datetime(1858, 11, 17)).days) + 3
    miss = [m for m in range(m0, m1 + 1) if m not in recs]
    if miss: sys.exit('ERP %s has no daily record for MJD %s' % (src, miss))
    with open(dst, 'w') as o:
        o.write('%% conveop                                 %s\n' % dt.datetime.utcnow().strftime('%d-%b-%y'))
        o.write('%%%% source: %s\n+pole&ut1\n' % os.path.basename(src))
        o.write('%% UT1 type = UT1\n%% Start&End%%Interval = %8d %8d   1.00\n' % (m0, m1))
        o.write('% Num. of Vars&Units =     5  0.1D+01  0.1D+01  0.1D+01  0.1D-02  0.1D-02\n')
        o.write('% Format = (f9.2,1x,2f10.6,f15.7,2f10.3,6(1x,a1))\n')
        o.write('%% MJD        XPOLE     YPOLE      UT1-TAI        DPSI     DEPSI    PRED_ID\n')
        for m in range(m0, m1 + 1):
            x, y, u = recs[m]
            o.write(' %8.2f %10.6f%10.6f %14.7f %9.3f %9.3f P P P\n' % (m, x, y, u, 0.0, 0.0))
        o.write('-pole&ut1\n')


def write_leap(dst):
    # TAI-UTC table in the GREAT format; the bundled file carries a non-existent 'predicted' leap second at 2021-01-01
    # (only the presence of the file matters to great_podflt, gtime has its own table, but keep it correct)
    rows = [(41317, 10), (41499, 11), (41683, 12), (42048, 13), (42413, 14), (42778, 15), (43144, 16), (43509, 17), (43874, 18),
            (44239, 19), (44786, 20), (45151, 21), (45516, 22), (46247, 23), (47161, 24), (47892, 25), (48257, 26), (48804, 27),
            (49169, 28), (49534, 29), (50083, 30), (50630, 31), (51179, 32), (53736, 33), (54832, 34), (56109, 35), (57204, 36),
            (57754, 37)]
    with open(dst, 'w') as o:
        o.write('+leap sec\n')
        for m, l in rows: o.write(' %5d  %2d\n' % (m, l))
        o.write('-leap sec\n')


def write_zero_bias(dst):
    # great_podflt needs a code-bias input (ALLBIAS); zero OSBs on the signals in use = no correction. IF clocks absorb
    # the common part; replace with a real SINEX-BIAS (e.g. CODE/GFZ OSB) for clocks consistent with the IGS convention.
    sig = {'G': 'C1C C1W C2W C2L C2S C2X C5Q C5X', 'E': 'C1C C1X C5Q C5X C7Q C7X C6C', 'C': 'C2I C6I C7I C1P C5P C7D',
           'R': 'C1C C1P C2C C2P', 'J': 'C1C C2L C5Q'}
    L = ['%=BIA 1.00 GRT 2026:001:00000 GRT 2000:001:00000 2099:365:00000 A 00000000', '+BIAS/SOLUTION',
         '*BIAS SVN_ PRN STATION__ OBS1 OBS2 BIAS_START____ BIAS_END______ UNIT __ESTIMATED_VALUE____ _STD_DEV___']
    for sysc, obs in sig.items():
        for n in range(1, 64 if sysc == 'C' else 37):
            for o in obs.split():
                L.append(' OSB  %s%03d %s%02d           %-4s      2000:001:00000 2099:365:00000 ns   %21.14E %11.5E '
                         % (sysc, n, sysc, n, o, 0.0, 0.1))
    L += ['-BIAS/SOLUTION', '%=ENDBIA']
    with open(dst, 'w') as f: f.write('\n'.join(L) + '\n')


def fetch_cod_osb(A, days):
    # CODE rapid OSB (absolute SINEX-BIAS, one file per day, ~1 day latency) from the IGS mirror. A day without its own
    # file (NRT) takes the newest earlier file within 7 days; its satellite records are re-dated to that day, because
    # great_podflt ignores (= zero) biases outside their validity window. Station records are dropped (not used).
    bdir = os.path.join(A.work, 'bias'); os.makedirs(bdir, exist_ok=True)
    def get(day):
        f = 'COD0OPSRAP_%s0000_01D_01D_OSB.BIA' % day.strftime('%Y%j'); p = os.path.join(bdir, f)
        if os.path.exists(p) and os.path.getsize(p) > 0: return p
        wk = (day - dt.datetime(1980, 1, 6)).days // 7
        for w in (wk, wk + 1):
            r = subprocess.run(['curl', '-sf', '-m', '300', '%s/%d/%s.gz' % (A.bias_mirror, w, f), '-o', p + '.gz'])
            if r.returncode == 0 and os.path.getsize(p + '.gz') > 0:
                subprocess.run(['gunzip', '-f', p + '.gz']); return p
            if os.path.exists(p + '.gz'): os.remove(p + '.gz')
        return None
    out = os.path.join(bdir, 'cod_osb_%s_%s.bia' % (days[0].strftime('%Y%j'), days[-1].strftime('%Y%j')))
    lines = []
    for day in days:
        src = None
        for back in range(0, 8):
            src = get(day - dt.timedelta(days=back))
            if src: break
        if not src: sys.exit('no CODE OSB file within 7 days before %s' % day.date())
        b = '%s:%03d:00000' % (day.year, day.timetuple().tm_yday)
        e1 = day + dt.timedelta(days=1); e = '%s:%03d:00000' % (e1.year, e1.timetuple().tm_yday)
        n = 0
        for ln in open(src, encoding='latin-1'):
            if ln.startswith(' OSB ') and ln[15:24].strip() == '':
                out_ln = (ln[:35] + b + ' ' + e + ln[64:].rstrip('\r\n')).ljust(103)
                lines.append(out_ln); n += 1
        print('bias %s: %s (%d records)' % (day.date(), os.path.basename(src), n))
    if not lines:
        sys.exit('no satellite OSB records found in downloaded CODE OSB files')
    with open(out, 'w') as f:
        f.write('%%=BIA 1.00 COD %s COD %s %s A %08d\n' % (lines[0][35:49], lines[0][35:49], lines[-1][50:64], len(lines)))
        f.write('+BIAS/SOLUTION\n' + '\n'.join(lines) + '\n-BIAS/SOLUTION\n%=ENDBIA\n')
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--stations', required=True); ap.add_argument('--beg', required=True); ap.add_argument('--end', required=True)
    ap.add_argument('--work', required=True); ap.add_argument('--tag', default=None)
    ap.add_argument('--rtcm-root', default=os.path.expanduser('~/gnss/rtcm'))
    ap.add_argument('--erp', required=True)
    ap.add_argument('--sys-dir', default=os.path.join(REPO, 'utils', 'gnss_sys_file'))
    ap.add_argument('--atx', default=None); ap.add_argument('--satinfo', default=None)
    ap.add_argument('--convbin', default=os.path.expanduser('~/.PRIDE_PPPAR_BIN/convbin'))
    ap.add_argument('--nav-streams', default='', help='extra stations whose streams relay ephemerides (comma list)')
    ap.add_argument('--sys', default='GPS', help='GPS, GAL, BDS, GLO or a space separated list')
    ap.add_argument('--max-stations', type=int, default=0); ap.add_argument('--min-epochs', type=int, default=240)
    ap.add_argument('--ref-clk', default=None, help='station ID (12 char) used as reference clock')
    ap.add_argument('--threads', type=int, default=8); ap.add_argument('-c', type=int, default=8)
    ap.add_argument('--force', action='store_true')
    ap.add_argument('--bias', default=None, help='SINEX-BIAS / Bernese DCB file (default: zero OSBs)')
    ap.add_argument('--bias-fetch', action='store_true', help='download CODE rapid OSB (COD0OPSRAP) for the arc days')
    ap.add_argument('--bias-mirror', default='ftp://igs.gnsswhu.cn/pub/gps/products', help='<mirror>/<GPS week>/<file>.gz')
    ap.add_argument('--leap', type=int, default=37, help='TAI-UTC leap seconds (default: 37)')
    ap.add_argument('--int', type=int, default=30, help='processing interval in seconds (default: 30)')
    ap.add_argument('--template', default=os.path.join(REPO, 'config', 'podflt_geodnet.xml.in'))
    ap.add_argument('--otl-models', default=os.environ.get('OTL_DIR'),
                    help='directory with FES2004S1.dat + Love_load_cm.dat: writes <work>/geodnet.blq with\n'
                         'scripts/make_blq.py (FES2004 CM frame); default $OTL_DIR')
    A = ap.parse_args()
    if A.int <= 0:
        sys.exit('--int must be positive')
    if not (shutil.which(A.convbin) or os.path.isfile(A.convbin)):
        sys.exit('convbin not found: %s (pass --convbin <path>)' % A.convbin)
    A.work = os.path.abspath(os.path.expanduser(A.work))
    if os.path.commonpath([A.work, REPO]) == REPO:
        sys.exit('--work must be outside the repository (it holds private station data)')
    beg, end = parse_t(A.beg), parse_t(A.end)
    if end <= beg:
        sys.exit('--end must be after --beg')
    os.makedirs(os.path.join(A.work, 'tmp'), exist_ok=True); os.chmod(A.work, 0o700)

    st = [(r['name'], (float(r['x']), float(r['y']), float(r['z']))) for r in csv.DictReader(open(A.stations))]
    if A.max_stations: st = st[:A.max_stations]
    amap = load_alias(os.path.join(A.work, 'alias_map.csv'), [n for n, _ in st])

    days = []
    d = dt.datetime(beg.year, beg.month, beg.day)
    while d < end: days.append(d); d += dt.timedelta(days=1)

    # RINEX observations (per station and day) and broadcast navigation (per day)
    jobs = [(n, amap[n], d) for n, _ in st for d in days]
    for d in days: os.makedirs(os.path.join(A.work, 'obs', d.strftime('%Y'), d.strftime('%j')), exist_ok=True)
    cnt, ok_days = {}, {}
    with cf.ThreadPoolExecutor(A.c) as ex:
        for (n, a, d), r in zip(jobs, ex.map(lambda j: to_rinex_obs(A, j[0], j[1], j[2], beg, end), jobs)):
            cnt[r] = cnt.get(r, 0) + 1
            if r in ('ok', 'skip', 'kept'): ok_days.setdefault(n, set()).add(d)
    print('RINEX obs:', cnt)
    for d in days: to_rinex_nav(A, d, beg, end)

    # stations used: RINEX for every day of the arc
    use = [(n, xyz) for n, xyz in st if len(ok_days.get(n, ())) == len(days)]
    if len(use) < 10: sys.exit('only %d stations with data for the whole arc' % len(use))
    print('stations with data on every day: %d of %d' % (len(use), len(st)))

    sites = [(amap[n], xyz) for n, xyz in use]
    write_sinex(os.path.join(A.work, 'sites.snx'), sites, beg, end)
    write_erp(A.erp, os.path.join(A.work, 'erp_great'), beg, end, leap=A.leap)
    write_leap(os.path.join(A.work, 'leap_seconds'))
    if A.bias_fetch:
        A.bias = fetch_cod_osb(A, days)
    elif not A.bias:
        A.bias = os.path.join(A.work, 'zero_osb.bia'); write_zero_bias(A.bias)

    blq = os.path.join(os.path.abspath(A.sys_dir), 'oceanload')
    if A.otl_models:
        blq = os.path.join(A.work, 'geodnet.blq')
        mid = (beg + (end - beg) / 2).strftime('%Y-%m-%d')
        subprocess.run([sys.executable, os.path.join(REPO, 'scripts', 'make_blq.py'), os.path.join(A.work, 'sites.snx'),
                        A.otl_models, blq, mid], check=True)
        os.chmod(blq, 0o600)
    else:
        print('warning: no --otl-models / $OTL_DIR: bundled oceanload (IGS sites only), no ocean loading for GEODNET sites')
    ref = amap[A.ref_clk] if A.ref_clk and A.ref_clk in amap else sites[0][0]
    tag = A.tag or '%s_%s' % (beg.strftime('%Y%m%dT%H'), end.strftime('%Y%m%dT%H'))
    run = os.path.join(A.work, 'run_' + tag)
    for sub in ('result', 'logger'): os.makedirs(os.path.join(run, sub), exist_ok=True)
    sysd = os.path.abspath(A.sys_dir)
    end_inc = max(beg, end - dt.timedelta(seconds=A.int))
    step_in = next(d for d in range(min(60, A.int), 0, -1) if A.int % d == 0)
    sub = {
        '@BEG@': beg.strftime('%Y-%m-%d %H:%M:%S'), '@END@': end_inc.strftime('%Y-%m-%d %H:%M:%S'), '@SYS@': A.sys,
        '@INT@': str(A.int), '@STEP_IN@': str(step_in), '@STEP_OUT@': str(A.int),
        '@SITES@': '\n'.join(' '.join(a.lower() for a, _ in sites[i:i + 10]) for i in range(0, len(sites), 10)),
        '@REFCLK@': ref.upper(), '@THREADS@': str(A.threads), '@RUN@': run, '@WORK@': A.work, '@SYSDIR@': sysd,
        '@ATX@': os.path.abspath(A.atx) if A.atx else os.path.join(sysd, 'igs20.atx'),
        '@SATINFO@': os.path.abspath(A.satinfo) if A.satinfo else os.path.join(sysd, 'satinfo'),
        '@SP3@': os.path.join(run, 'result', 'orbit.sp3'), '@BIAS@': os.path.abspath(A.bias), '@BLQ@': blq,
    }
    xml = open(A.template).read()
    for k, v in sub.items(): xml = xml.replace(k, v)
    with open(os.path.join(run, 'podflt.xml'), 'w') as f: f.write(xml)
    shutil.rmtree(os.path.join(A.work, 'tmp'), ignore_errors=True)
    print('run directory:', run)


if __name__ == '__main__':
    main()
