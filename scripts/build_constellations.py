#!/usr/bin/env python3
"""
Builds constellations.json — pre-computed line segments for 10 well-known
constellation figures, cross-referenced against HYG so the endpoint
positions match the star catalog exactly.

To add or remove figures: edit CONSTELLATIONS below. Each line is a pair
(bayer_letter, constellation_abbrev) — e.g. ('Alp', 'UMa') for Alpha
Ursae Majoris. We look that up in HYG's `bf` column, prefer the primary
component for double stars (Alp1Cen, not Alp2Cen), and drop any star
beyond R_MAX_PC.

Input:  data/hyg_v42.csv  (see scripts/build_stars.py for download URL)
Output: constellations.json — { r_max, figures: [{name, segments}, …] }
"""
import csv, json, math, re, sys
from pathlib import Path

HERE     = Path(__file__).resolve().parent
ROOT     = HERE.parent
DATA_DIR = ROOT / 'data'
HYG_PATH = DATA_DIR / 'hyg_v42.csv'
OUT_PATH = ROOT / 'constellations.json'

R_MAX_PC = 200.0    # drop a figure's segment if either endpoint is farther than this

# Each figure: list of (bayer_a, con_a, bayer_b, con_b) line endpoints.
# Bayer is the 3-letter Greek abbrev; constellation abbrev is the IAU 3-letter
# code (UMa, UMi, Cas, Cru, …).
CONSTELLATIONS = [
    {'name': 'URSA MAJOR', 'lines': [
        ('Eta','UMa','Zet','UMa'), ('Zet','UMa','Eps','UMa'),
        ('Eps','UMa','Del','UMa'), ('Del','UMa','Gam','UMa'),
        ('Gam','UMa','Bet','UMa'), ('Bet','UMa','Alp','UMa'),
        ('Alp','UMa','Del','UMa'),                            # bowl-close diagonal
    ]},
    {'name': 'URSA MINOR', 'lines': [
        ('Alp','UMi','Del','UMi'), ('Del','UMi','Eps','UMi'),
        ('Eps','UMi','Zet','UMi'), ('Zet','UMi','Eta','UMi'),
        ('Eta','UMi','Bet','UMi'), ('Bet','UMi','Gam','UMi'),
        ('Gam','UMi','Zet','UMi'),
    ]},
    {'name': 'CASSIOPEIA', 'lines': [
        ('Eps','Cas','Del','Cas'), ('Del','Cas','Gam','Cas'),
        ('Gam','Cas','Alp','Cas'), ('Alp','Cas','Bet','Cas'),
    ]},
    {'name': 'BOÖTES', 'lines': [
        ('Alp','Boo','Eps','Boo'), ('Eps','Boo','Gam','Boo'),
        ('Gam','Boo','Bet','Boo'), ('Bet','Boo','Del','Boo'),
        ('Del','Boo','Alp','Boo'), ('Alp','Boo','Eta','Boo'),  # Muphrid spur
    ]},
    {'name': 'LEO', 'lines': [
        # Sickle
        ('Alp','Leo','Gam','Leo'), ('Gam','Leo','Zet','Leo'),
        ('Zet','Leo','Mu','Leo'),  ('Mu','Leo','Eps','Leo'),
        # Body triangle
        ('Alp','Leo','The','Leo'), ('The','Leo','Bet','Leo'),
        ('The','Leo','Del','Leo'), ('Del','Leo','Gam','Leo'),
    ]},
    {'name': 'PEGASUS', 'lines': [
        ('Alp','Peg','Bet','Peg'), ('Bet','Peg','Alp','And'),  # shared corner w/ Andromeda
        ('Alp','And','Gam','Peg'), ('Gam','Peg','Alp','Peg'),
    ]},
    {'name': 'AQUILA', 'lines': [
        ('Gam','Aql','Alp','Aql'), ('Alp','Aql','Bet','Aql'),
        ('Alp','Aql','Zet','Aql'), ('Alp','Aql','Del','Aql'),
        ('Del','Aql','The','Aql'), ('Del','Aql','Lam','Aql'),
    ]},
    {'name': 'HERCULES', 'lines': [
        ('Eta','Her','Zet','Her'), ('Zet','Her','Eps','Her'),
        ('Eps','Her','Pi','Her'),  ('Pi','Her','Eta','Her'),
    ]},
    {'name': 'CRUX', 'lines': [
        ('Alp','Cru','Gam','Cru'), ('Bet','Cru','Del','Cru'),
    ]},
    {'name': 'CENTAURUS', 'lines': [
        ('Alp','Cen','Bet','Cen'), ('Bet','Cen','Eps','Cen'),
        ('Eps','Cen','Zet','Cen'), ('Zet','Cen','Gam','Cen'),
    ]},
]

# HYG packs Flamsteed + Bayer + constellation into the `bf` column with format
# quirks: leading Flamsteed digits, optional super-script digit, double-space
# alignment for 2-letter Bayer (e.g. "34Mu  UMa", "41Gam1Leo", "Alp1Cen").
# Returns (bayer, super, con) or None for pure-Flamsteed/unparseable entries.
def parse_bf(bf):
    m = re.match(r'^(\d*)(.*)$', bf)
    if not m: return None
    rest = m.group(2)
    # Constellation abbrev: first letter upper, second either case (handles UMa/CMa/TrA), third lower.
    m2 = re.match(r'^(.+?)\s*([A-Z][A-Za-z][a-z])$', rest)
    if not m2: return None
    bayer_part = m2.group(1).strip()
    con = m2.group(2)
    if not bayer_part: return None
    super_ = ''
    if bayer_part[-1].isdigit():
        super_ = bayer_part[-1]
        bayer_part = bayer_part[:-1]
    return (bayer_part, super_, con)

def main():
    if not HYG_PATH.exists():
        print(f'ERROR: {HYG_PATH} not found. See data/README.md.', file=sys.stderr)
        sys.exit(1)

    with HYG_PATH.open() as f:
        rows = list(csv.DictReader(f))

    # (bayer, con) -> row; prefer primary component (no superscript or '1').
    idx, super_idx = {}, {}
    def rank(s): return 0 if s in ('', '1') else (1 if s == '2' else 2)
    for r in rows:
        bf = (r.get('bf') or '').strip()
        if not bf: continue
        p = parse_bf(bf)
        if not p: continue
        bayer, sup, con = p
        key = (bayer, con)
        if key not in idx or rank(sup) < rank(super_idx[key]):
            idx[key] = r; super_idx[key] = sup

    result, total_segs = [], 0
    missing, oor = [], []
    for c in CONSTELLATIONS:
        segs = []
        for ba, ca, bb, cb in c['lines']:
            ra = idx.get((ba, ca)); rb = idx.get((bb, cb))
            if ra is None or rb is None:
                missing.append((c['name'], ba, ca, bb, cb)); continue
            da, db = float(ra['dist']), float(rb['dist'])
            if da > R_MAX_PC or db > R_MAX_PC:
                oor.append((c['name'], f"{ba}{ca}" if da > R_MAX_PC else f"{bb}{cb}",
                            da if da > R_MAX_PC else db)); continue
            segs.append([
                round(float(ra['x']),3), round(float(ra['y']),3), round(float(ra['z']),3),
                round(float(rb['x']),3), round(float(rb['y']),3), round(float(rb['z']),3),
            ])
        total_segs += len(segs)
        result.append({'name': c['name'], 'segments': segs})
        print(f'  {c["name"]:14s} {len(segs):2d}/{len(c["lines"]):2d} lines')

    if missing:
        print('\nMISSING stars:')
        for m in missing: print(' ', m)
    if oor:
        print('\nOUT OF RANGE (dropped):')
        for o in oor: print(' ', o)

    OUT_PATH.write_text(
        json.dumps({'r_max': R_MAX_PC, 'figures': result},
                   separators=(',', ':'), ensure_ascii=False),
        encoding='utf-8',
    )
    print(f'\nWrote {OUT_PATH} ({OUT_PATH.stat().st_size} bytes), {total_segs} segments.')

if __name__ == '__main__':
    main()
