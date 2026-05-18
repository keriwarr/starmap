#!/usr/bin/env python3
"""
Builds stars.json for the starmap from a HYG + GCNS hybrid catalog.

Why hybrid:
  - HYG carries the bright/famous stars (Sirius, Vega, Polaris, …) plus
    proper names + Bayer/Flamsteed designations. Gaia EDR3 saturates above
    G ≈ 3, so those navigational beacons are missing from GCNS entirely.
  - GCNS fills in the volume-complete dim end out to 100 pc — lots of
    M dwarfs Gaia DR3 measures cleanly that HYG doesn't list.

Cross-catalog dedup is by sky direction, not 3-D distance. Same physical
star can sit several pc apart in the two catalogs because of parallax
disagreement, but the sky direction agrees to <1″ after proper-motion
correction. We bring HYG positions to J2016 (Gaia EDR3 epoch) using HYG's
pmra/pmdec (mas/yr; pmra is already μα* = μα·cos δ), then drop any GCNS
row whose unit-sphere position is within SKY_TOL_ARCSEC of an HYG entry.

Within-catalog dedup is intentionally avoided. HYG explicitly lists
Sirius A + B, α Cen A + B, etc. as separate rows — that's the catalog's
design, not noise.

Inputs (download once, put under data/):
  data/hyg_v42.csv   — HYG v4.2 main table (uncompressed)
      https://codeberg.org/astronexus/hyg/media/branch/main/data/hyg/CURRENT/hyg_v42.csv.gz
  data/gcns_table1c.dat.gz  — GCNS main table (Smart et al. 2020)
      https://cdsarc.u-strasbg.fr/ftp/J/A+A/649/A6/table1c.dat.gz

Output:
  stars.json — { meta: {r_1000, r_max, count}, stars: [[name, x, y, z, mag, spect], …] }
  All positions are heliocentric equatorial Cartesian in parsecs (HYG
  convention). `mag` is V-band for HYG entries, Gaia G for GCNS entries
  (close enough for visualization). `spect` is one of OBAFGKMLTYWC or '?'.
"""
import csv, gzip, json, math, os, re, sys
from collections import defaultdict
from pathlib import Path

HERE      = Path(__file__).resolve().parent
ROOT      = HERE.parent
DATA_DIR  = ROOT / 'data'
HYG_PATH  = DATA_DIR / 'hyg_v42.csv'
GCNS_PATH = DATA_DIR / 'gcns_table1c.dat.gz'
OUT_PATH  = ROOT / 'stars.json'

# --- knobs ------------------------------------------------------------------
HYG_R_MAX_PC    = 200.0     # keep HYG out to here (constellation backbones live here)
GCNS_R_MAX_PC   = 100.0     # GCNS's native horizon
GCNS_MAG_LIMIT  = 14.0      # Gaia G; bump up for more dim fillers (at file-size cost)
SKY_TOL_ARCSEC  = 2.0       # cross-catalog dedup tolerance after PM correction
EPOCH_DT_YEARS  = 16.0      # J2000 (HYG) -> J2016 (Gaia EDR3)
SKY_TOL_RAD     = SKY_TOL_ARCSEC / 3600.0 * math.pi / 180.0
SKY_CELL        = SKY_TOL_RAD  # cell size = tolerance; ±1 cell search covers it

# --- HYG --------------------------------------------------------------------
_PROPER_BLOCKLIST = re.compile(r'^(HIP|HD|Gl|GJ|ID|TYC|2MASS)\s', re.I)
_BAYER_FLAM_SUFFIX = re.compile(r'^[A-Z][a-z]{2}$')

def is_proper(name):
    """Heuristic: a name is 'proper' if it isn't a catalog ID, doesn't contain
    digits (no Flamsteed), and doesn't end in a 3-letter capitalized
    constellation abbreviation (no Bayer designations like 'Alp Cen')."""
    if not name: return False
    if name == 'Sol': return True
    if _PROPER_BLOCKLIST.match(name): return False
    if any(c.isdigit() for c in name): return False
    parts = name.split()
    if len(parts) >= 2 and _BAYER_FLAM_SUFFIX.match(parts[-1]):
        return False
    return True

def spect_class_hyg(s):
    if not s: return '?'
    m = re.match(r'(?:sd|D)?([OBAFGKMLTYWC])', s)
    return m.group(1) if m else '?'

def hyg_pm_corrected_uvec(r):
    """HYG row → unit-sphere position at J2016. Returns None if unparseable."""
    try:
        ra_h = float(r['ra']); dec_deg = float(r['dec'])
    except (TypeError, ValueError):
        return None
    try: pmra = float(r['pmra'] or 0)
    except ValueError: pmra = 0
    try: pmdec = float(r['pmdec'] or 0)
    except ValueError: pmdec = 0
    if abs(pmra)  > 1e5: pmra  = 0   # HYG sometimes uses sentinel-ish bad values
    if abs(pmdec) > 1e5: pmdec = 0
    ra_deg = ra_h * 15.0
    cosd = math.cos(math.radians(dec_deg))
    if abs(cosd) < 1e-6: cosd = 1e-6 if cosd >= 0 else -1e-6
    # HYG pmra is μα* (mas/yr, already × cos δ); recover dra/dt by dividing.
    ra_deg += (pmra / cosd) * EPOCH_DT_YEARS / 3.6e6
    dec_deg += pmdec * EPOCH_DT_YEARS / 3.6e6
    rar = math.radians(ra_deg); decr = math.radians(dec_deg)
    cd = math.cos(decr)
    return (cd * math.cos(rar), cd * math.sin(rar), math.sin(decr))

def load_hyg_and_sky_index(r_max):
    """Return (rows, sky_grid). rows are output-ready [name, x, y, z, mag, spect]."""
    rows = []
    grid = defaultdict(list)
    with HYG_PATH.open() as f:
        for r in csv.DictReader(f):
            try:
                d = float(r['dist'])
                if d > r_max: continue
                x = round(float(r['x']), 3)
                y = round(float(r['y']), 3)
                z = round(float(r['z']), 3)
                mag = float(r['mag']) if r['mag'] else None
            except (TypeError, ValueError):
                continue
            n = r['proper'] or ''
            n = n if is_proper(n) else ''
            rows.append([n, x, y, z,
                         round(mag, 2) if mag is not None else None,
                         spect_class_hyg(r['spect'])])
            u = hyg_pm_corrected_uvec(r)
            if u is None: continue
            ux, uy, uz = u
            grid[(int(ux/SKY_CELL), int(uy/SKY_CELL), int(uz/SKY_CELL))].append((ux, uy, uz))
    return rows, grid

# --- GCNS -------------------------------------------------------------------
def spect_from_bp_rp(bp_rp):
    """Approx MK class from Gaia BP-RP color. Boundaries match the renderer's palette."""
    if bp_rp is None: return '?'
    if bp_rp < -0.30: return 'O'
    if bp_rp <  0.00: return 'B'
    if bp_rp <  0.50: return 'A'
    if bp_rp <  0.75: return 'F'
    if bp_rp <  1.00: return 'G'
    if bp_rp <  1.50: return 'K'
    return 'M'

def radec_to_xyz(ra_deg, dec_deg, dist_pc):
    ra = math.radians(ra_deg); dec = math.radians(dec_deg)
    cd = math.cos(dec)
    return (dist_pc * cd * math.cos(ra),
            dist_pc * cd * math.sin(ra),
            dist_pc * math.sin(dec))

def radec_to_uvec(ra_deg, dec_deg):
    ra = math.radians(ra_deg); dec = math.radians(dec_deg)
    cd = math.cos(dec)
    return (cd * math.cos(ra), cd * math.sin(ra), math.sin(dec))

# GCNS table1c.dat column byte ranges (1-indexed, inclusive) per the ReadMe at
# https://cdsarc.u-strasbg.fr/ftp/J/A+A/649/A6/ReadMe . Fixed-width format.
def parse_gcns_row(line):
    def fld(a, b): return line[a-1:b].strip()
    try:
        ra_s    = fld(23, 36)
        dec_s   = fld(46, 59)
        plx_s   = fld(69, 77)
        g_s     = fld(123, 130)
        bp_s    = fld(142, 149)
        rp_s    = fld(161, 168)
        prob_s  = fld(240, 244)
        dist50  = fld(278, 289)  # kpc
        if not ra_s or not dec_s or not g_s: return None
        ra, dec, gmag = float(ra_s), float(dec_s), float(g_s)
        # Prefer Bayesian Dist50 (kpc → pc); fall back to 1/parallax.
        if dist50:
            dist_pc = float(dist50) * 1000.0
        elif plx_s:
            plx = float(plx_s)
            if plx <= 0: return None
            dist_pc = 1000.0 / plx
        else:
            return None
        if dist_pc <= 0 or dist_pc > GCNS_R_MAX_PC: return None
        bp_rp = None
        if bp_s and rp_s:
            try: bp_rp = float(bp_s) - float(rp_s)
            except ValueError: pass
        prob = float(prob_s) if prob_s else 0.0
        return (ra, dec, dist_pc, gmag, bp_rp, prob)
    except (ValueError, IndexError):
        return None

def has_sky_match(grid, ux, uy, uz):
    cx, cy, cz = int(ux/SKY_CELL), int(uy/SKY_CELL), int(uz/SKY_CELL)
    t2 = SKY_TOL_RAD * SKY_TOL_RAD
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                for hx, hy, hz in grid.get((cx+dx, cy+dy, cz+dz), ()):
                    if (hx-ux)**2 + (hy-uy)**2 + (hz-uz)**2 <= t2:
                        return True
    return False

# --- main -------------------------------------------------------------------
def main():
    for p in (HYG_PATH, GCNS_PATH):
        if not p.exists():
            print(f'ERROR: required file missing: {p}', file=sys.stderr)
            print('See data/README.md for download instructions.', file=sys.stderr)
            sys.exit(1)

    print(f'HYG: loading + PM-correcting sky positions within {HYG_R_MAX_PC} pc...')
    hyg, sky_grid = load_hyg_and_sky_index(HYG_R_MAX_PC)
    named = sum(1 for r in hyg if r[0])
    print(f'  {len(hyg)} HYG stars ({named} proper-named).')

    print(f'GCNS: streaming table1c, sky dedup at {SKY_TOL_ARCSEC}", mag G < {GCNS_MAG_LIMIT}...')
    gcns_only = []
    kept = matched = rejected = 0
    with gzip.open(GCNS_PATH, 'rt', encoding='latin-1') as f:
        for line in f:
            p = parse_gcns_row(line)
            if p is None: rejected += 1; continue
            ra, dec, dist_pc, gmag, bp_rp, prob = p
            if prob < 0.5 or gmag > GCNS_MAG_LIMIT:
                rejected += 1; continue
            ux, uy, uz = radec_to_uvec(ra, dec)
            if has_sky_match(sky_grid, ux, uy, uz):
                matched += 1; continue
            x, y, z = radec_to_xyz(ra, dec, dist_pc)
            gcns_only.append(['', round(x, 3), round(y, 3), round(z, 3),
                              round(gmag, 2), spect_from_bp_rp(bp_rp)])
            kept += 1
    print(f'  GCNS kept: {kept}, dropped-as-HYG-dup: {matched}, otherwise rejected: {rejected}')

    combined = hyg + gcns_only
    def dist(row): return math.sqrt(row[1]**2 + row[2]**2 + row[3]**2)
    combined.sort(key=dist)

    r_1000 = dist(combined[1000]) if len(combined) > 1000 else dist(combined[-1])
    r_max  = dist(combined[-1])

    payload = {
        'meta': {'r_1000': round(r_1000, 3),
                 'r_max':  round(r_max, 3),
                 'count':  len(combined)},
        'stars': combined,
    }
    OUT_PATH.write_text(json.dumps(payload, separators=(',', ':'), ensure_ascii=False),
                        encoding='utf-8')

    final_named = sum(1 for s in combined if s[0])
    sz_kb = OUT_PATH.stat().st_size / 1024
    print(f'\nWrote {OUT_PATH}')
    print(f'  {len(combined)} stars ({final_named} named)')
    print(f'  R_1000 = {r_1000:.3f} pc, R_max = {r_max:.3f} pc')
    print(f'  raw: {sz_kb:.1f} KB')

if __name__ == '__main__':
    main()
