# Build scripts

Two Python scripts regenerate the catalog files the renderer loads:

- `build_stars.py` → `stars.json` (the main star catalog, HYG + GCNS hybrid)
- `build_constellations.py` → `constellations.json` (10 figure-line graphs cross-referenced against HYG)

## Setup

Both scripts read raw catalogs from `data/`. Create the directory and pull
the two source files (kept out of git — large + redistributable):

```bash
mkdir -p data

# HYG v4.2 — Hipparcos + Yale Bright Star + Gliese merge (~34 MB after gunzip)
curl -L -o data/hyg_v42.csv.gz \
  https://codeberg.org/astronexus/hyg/media/branch/main/data/hyg/CURRENT/hyg_v42.csv.gz
gunzip data/hyg_v42.csv.gz

# GCNS — Gaia Catalogue of Nearby Stars, Smart et al. 2021 (~75 MB, leave compressed)
curl -kL -o data/gcns_table1c.dat.gz \
  https://cdsarc.u-strasbg.fr/ftp/J/A+A/649/A6/table1c.dat.gz
```

## Run

```bash
python3 scripts/build_stars.py
python3 scripts/build_constellations.py
```

No third-party Python dependencies — stdlib only.

## Knobs (top of `build_stars.py`)

| Constant | Current | Effect |
|---|---|---|
| `HYG_R_MAX_PC` | 200 | Outer radius for HYG. Constellation figures live out here. |
| `GCNS_R_MAX_PC` | 100 | GCNS's native completeness horizon. |
| `GCNS_MAG_LIMIT` | 14.0 | Gaia G cutoff. Higher = more dim fillers, bigger JSON. |
| `SKY_TOL_ARCSEC` | 2.0 | Cross-catalog dedup tolerance. Real binaries with wider separation stay. |

## Pipeline notes

- **Cross-catalog dedup is by sky direction, not 3-D distance.** Same physical
  star can sit several pc apart between the two catalogs because of parallax
  disagreement, but sky directions agree to <1″ once HYG is proper-motion-
  corrected to J2016 (Gaia EDR3 epoch). 3-D dedup misses most matches and
  produces a characteristic "rays from Sol" pattern as duplicates render at
  the same direction but different distances.
- **HYG is preserved verbatim** (no intra-HYG dedup). Sirius A + B, α Cen A
  + B, etc. are real catalog entries.
- **No intra-GCNS dedup.** Gaia legitimately resolves some binaries into
  separate source IDs; both stay.
- HYG `pmra` is `μα* = μα·cos δ` in mas/yr (Hipparcos convention) — divide by
  `cos δ` to recover `dra/dt`.
- HYG `ra` is in HOURS (not degrees) — multiply by 15.
- HYG `bf` packs Flamsteed + Bayer + IAU abbrev (e.g. `34Mu  UMa`,
  `41Gam1Leo`, `Alp1Cen`). The parser in `build_constellations.py` handles
  the quirks.

## Output schemas

`stars.json`:
```jsonc
{
  "meta": { "r_1000": 12.907, "r_max": 200.0, "count": 146392 },
  "stars": [
    ["Sol", 0.0, 0.0, 0.0, -26.7, "G"],
    ["Proxima Centauri", -0.472, -0.361, -1.151, 11.01, "M"],
    ["", -0.495, -0.414, -1.157, 1.35, "K"],   // unnamed: 6-element row, name=""
    ...
  ]
}
```

`constellations.json`:
```jsonc
{
  "r_max": 200.0,
  "figures": [
    { "name": "URSA MAJOR",
      "segments": [[x1,y1,z1,x2,y2,z2], ...] },
    ...
  ]
}
```

`x, y, z` are heliocentric equatorial Cartesian in parsecs (HYG convention,
X toward vernal equinox). `mag` is V-band for HYG / Gaia G for GCNS rows.
`spect` is the MK class letter (O B A F G K M …) or `?`.
