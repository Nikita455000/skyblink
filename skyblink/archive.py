"""Data access: SPHEREx spectral images and Gaia DR3 from NASA/IPAC IRSA, names from CDS Sesame."""
import concurrent.futures
import csv
import hashlib
import io
import json
import logging
import math
import os
import re
import tempfile
import time
import urllib.parse
import warnings

import numpy as np
import requests
from astropy.io import fits
from astropy.wcs import WCS
from astroquery.ipac.irsa import Irsa

logging.getLogger("astropy").setLevel(logging.ERROR)
warnings.filterwarnings("ignore", module="astropy")

IRSA = "https://irsa.ipac.caltech.edu/"
SESAME = "https://cds.unistra.fr/cgi-bin/nph-sesame/-oxp/SNV?"
CACHE_DIR = os.environ.get("SKYBLINK_CACHE") or os.path.join(tempfile.gettempdir(), "skyblink-cache")
BLOCK = 2880

# FLAGS bits (names from the FLAGS header of SPHEREx images). OUTLIER is never masked on purpose:
# a moving object is an "outlier" against the static sky.
# Always unusable: SUR_ERROR, NONFUNC, MISSING_DATA, HOT, COLD.
BAD_MASK = sum(1 << b for b in (2, 6, 9, 10, 11))
# Fine to look at but not to measure: TRANSIENT, OVERFLOW, NONLINEAR (saturated stars), STREAK, BLOOM, SNOWBALL.
SOFT_MASK = sum(1 << b for b in (0, 1, 15, 25, 26, 27))
CACHE_VERSION = "v2"

# Wavelength coverage of the six SPHEREx detector arrays, in microns.
DETECTORS = {1: (0.75, 1.09), 2: (1.10, 1.62), 3: (1.63, 2.41), 4: (2.42, 3.82), 5: (3.83, 4.41), 6: (4.42, 5.00)}

_session = requests.Session()
_session.headers["User-Agent"] = "SkyBlink/1.0 (NASA Space Apps 2026, team DeepSpace)"


def _cache_path(name):
    os.makedirs(CACHE_DIR, exist_ok=True)
    return os.path.join(CACHE_DIR, name)


# ---------------------------------------------------------------- names


# The challenge's test targets, so a demo never waits on a name resolver.
# J2000 positions and proper motions (mas/yr) as returned by CDS Sesame / SIMBAD.
KNOWN = {
    "barnard's star": ("Barnard's Star", 269.4520769586187, 4.69336496657667, -801.551, 10362.394),
    "eps ind b": ("Epsilon Indi B", 331.04414658852124, -56.78284419511611, 3981.977, -2466.832),
    "wise j085510.83-071442.5": ("WISE 0855−0714", 133.81890825, -7.247048, -8123.7, 673.2),
    "rw aur": ("RW Aurigae", 76.95652603458, 30.40143817389, 7.23, -20.33),
    "r dra": ("R Draconis", 248.16761083295, 66.75497941722, -20.643, 6.149),
}
ALIASES = {"barnard": "barnard's star", "barnards star": "barnard's star", "epsilon indi b": "eps ind b",
           "wise 0855": "wise j085510.83-071442.5", "wise 0855-0714": "wise j085510.83-071442.5",
           "rw aurigae": "rw aur", "r draconis": "r dra"}
SESAME_MIRRORS = [SESAME, "https://vizier.cfa.harvard.edu/viz-bin/nph-sesame/-oxp/SNV?"]


def resolve_name(name):
    """Resolve an object name: built-in challenge targets first, then CDS Sesame (with a mirror).

    Returns the J2000 position and proper motion (mas/yr), or None if the name is unknown.
    """
    key = " ".join(name.lower().split())
    key = ALIASES.get(key, key)
    if key in KNOWN:
        label, ra, dec, pmra, pmdec = KNOWN[key]
        return {"name": label, "ra": ra, "dec": dec, "pmra": pmra, "pmdec": pmdec, "epoch": 2000.0}
    text = None
    for base in SESAME_MIRRORS:
        try:
            r = _session.get(base + urllib.parse.quote(name.strip()), timeout=15)
            r.raise_for_status()
            text = r.text
            break
        except requests.RequestException:
            continue
    if text is None:
        raise IOError("name resolver unavailable")

    def tag(t):
        m = re.search(f"<{t}>([^<]+)</{t}>", text)
        return m.group(1).strip() if m else None

    if tag("jradeg") is None:
        return None
    return {
        "name": (tag("oname") or name).replace("NAME ", "").replace("V* ", "").replace("* ", ""),
        "ra": float(tag("jradeg")),
        "dec": float(tag("jdedeg")),
        "pmra": float(tag("pmRA") or 0.0),
        "pmdec": float(tag("pmDE") or 0.0),
        "epoch": 2000.0,
    }


# ---------------------------------------------------------------- SPHEREx image search


def _tap(query, tries=3):
    for k in range(tries):
        try:
            return Irsa.query_tap(query).to_table()
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(2 * (k + 1))


def query_images(ra, dec, max_age_hours=12):
    """All SPHEREx spectral images that cover (ra, dec), oldest first. Cached on disk for a few hours."""
    path = _cache_path(f"tap_{ra:.4f}_{dec:+.4f}.json")
    if os.path.exists(path) and time.time() - os.path.getmtime(path) < max_age_hours * 3600:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    t = _tap(f"""
        SELECT a.uri, p.time_bounds_lower, p.energy_bandpassname
        FROM spherex.artifact a
        JOIN spherex.plane p ON a.planeid = p.planeid
        WHERE 1 = CONTAINS(POINT('ICRS', {ra:.6f}, {dec:.6f}), p.poly)
        ORDER BY p.time_bounds_lower""")
    items = [{"uri": str(u), "mjd": float(m), "detector": int(str(b)[-1])}
             for u, m, b in zip(t["uri"], t["time_bounds_lower"], t["energy_bandpassname"])]
    with open(path, "w", encoding="utf-8") as f:
        json.dump(items, f)
    return items


# ---------------------------------------------------------------- streamed FITS cutouts


def _read_exact(raw, n):
    out = bytearray()
    while len(out) < n:
        chunk = raw.read(n - len(out))
        if not chunk:
            raise IOError("truncated FITS stream")
        out += chunk
    return bytes(out)


def _read_header(raw):
    """Read one FITS header (whole 2880-byte blocks up to END). None at a clean end of stream."""
    first = raw.read(BLOCK)
    if not first:
        return None
    buf = first if len(first) == BLOCK else first + _read_exact(raw, BLOCK - len(first))
    while True:
        blk = buf[-BLOCK:]
        if any(blk[i:i + 8] == b"END     " for i in range(0, BLOCK, 80)):
            return buf
        buf += _read_exact(raw, BLOCK)


def _data_bytes(header):
    naxis = header.get("NAXIS", 0)
    if naxis == 0:
        return 0
    n = abs(header["BITPIX"]) // 8 * math.prod(header[f"NAXIS{i}"] for i in range(1, naxis + 1))
    n += header.get("PCOUNT", 0)
    return int(math.ceil(n / BLOCK) * BLOCK)


def _download(url, want=None, tries=8):
    """Stream a FITS file, stopping as soon as the HDUs named in `want` have arrived.

    SPHEREx cutouts end with a large PSF table we never need, so stopping early saves ~4 MB per image.
    IRSA answers "Server busy" under load, so failures are retried with a growing pause.
    """
    for k in range(tries):
        try:
            with _session.get(url, stream=True, timeout=120) as r:
                if r.status_code != 200:
                    raise IOError(f"HTTP {r.status_code}")
                chunks, got, primary = [], set(), True
                while True:
                    hb = _read_header(r.raw)
                    if hb is None:
                        break
                    if primary and not hb.startswith(b"SIMPLE"):
                        raise IOError("not a FITS file (server busy)")
                    h = fits.Header.fromstring(hb.decode("ascii", "replace"))
                    chunks.append(hb + _read_exact(r.raw, _data_bytes(h)))
                    if not primary:
                        got.add(h.get("EXTNAME"))
                    primary = False
                    if want and set(want) <= got:
                        break
                if want and not set(want) <= got:
                    raise IOError("missing extensions")
                return b"".join(chunks)
        except (IOError, requests.RequestException):
            time.sleep(min(1.5 * (k + 1), 8))
    return None


# Keywords needed for the sky (TAN-SIP) and wavelength (WAVE-TAB, key "W") coordinate systems.
_WCS_KEY = re.compile(r"^(NAXIS\d?|WCSAXES\w?|WCSNAME\w?|CTYPE\d\w?|CRVAL\d\w?|CRPIX\d\w?|CDELT\d\w?|CUNIT\d\w?|"
                      r"CNAME\d\w?|PC\d_\d\w?|CD\d_\d\w?|PV\d_\d+\w?|PS\d_\d+\w?|[AB]P?_(ORDER|\d+_\d+)|"
                      r"RADESYS\w?|EQUINOX\w?|LONPOLE\w?|LATPOLE\w?)$")


def _wcs_header(header):
    """Keep only coordinate keywords: building a WCS from the full 270-card header is ~20x slower."""
    return fits.Header([(c.keyword, c.value) for c in header.cards if _WCS_KEY.match(c.keyword)]).tostring()


def fetch_cutout(item, ra, dec, size_deg):
    """One SPHEREx cutout as {data, soft, header, mjd, detector}.

    data is IMAGE minus ZODI with unusable pixels set to NaN; soft marks pixels to skip when measuring fluxes.
    """
    key = hashlib.sha1(f"{CACHE_VERSION}|{item['uri']}|{ra:.5f}|{dec:.5f}|{size_deg:.4f}".encode()).hexdigest()
    path = _cache_path(key + ".npz")
    meta = {k: item[k] for k in ("uri", "mjd", "detector")}
    if os.path.exists(path):
        z = np.load(path)
        return {"data": z["data"], "soft": z["soft"], "header": str(z["header"]), **meta}
    raw = _download(f"{IRSA}{item['uri']}?center={ra:.6f},{dec:.6f}&size={size_deg:.4f}",
                    want=("IMAGE", "FLAGS", "ZODI"))
    if raw is None:
        return None
    with fits.open(io.BytesIO(raw)) as hl:
        data = hl["IMAGE"].data.astype(np.float32) - hl["ZODI"].data.astype(np.float32)
        flags = hl["FLAGS"].data
        data[(flags & BAD_MASK) != 0] = np.nan
        soft = (flags & SOFT_MASK) != 0
        header = _wcs_header(hl["IMAGE"].header)
    np.savez_compressed(path, data=data, soft=soft, header=np.array(header))
    return {"data": data, "soft": soft, "header": header, **meta}


def fetch_cutouts(items, ra, dec, size_deg, progress=None, workers=12):
    """Download many cutouts in parallel. `progress(done, total)` is called as they arrive."""
    out = [None] * len(items)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(fetch_cutout, it, ra, dec, size_deg): i for i, it in enumerate(items)}
        for done, f in enumerate(concurrent.futures.as_completed(futures), 1):
            try:
                out[futures[f]] = f.result()
            except Exception:
                out[futures[f]] = None
            if progress:
                progress(done, len(items))
    return [c for c in out if c is not None]


# ---------------------------------------------------------------- wavelength lookup


def wave_table(item, ra, dec):
    """The WCS-WAVE lookup table for this image's detector (identical for all images of a detector)."""
    version = (re.search(r"l2b-v\d+", item["uri"]) or [""])[0]
    name = f"wcswave_D{item['detector']}_{version}.fits"
    bundled = os.path.join(os.path.dirname(__file__), "data", name)
    path = bundled if os.path.exists(bundled) else _cache_path(name)
    if not os.path.exists(path):
        raw = _download(f"{IRSA}{item['uri']}?center={ra:.6f},{dec:.6f}&size=0.01", want=("WCS-WAVE",))
        if raw is None:
            return None
        with fits.open(io.BytesIO(raw)) as hl:
            fits.HDUList([fits.PrimaryHDU(), hl["WCS-WAVE"].copy()]).writeto(path, overwrite=True)
    with fits.open(path) as hl:
        return hl["WCS-WAVE"].copy()


def wavelength(header, table, x, y):
    """Central wavelength (micron) of cutout pixel (x, y), zero-based.

    Evaluates the WAVE-TAB lookup (FITS Paper III) directly: building an astropy WCS for it
    costs ~30 ms per image, this costs microseconds. `wavelength_astropy` is the reference.
    """
    hdr = fits.Header.fromstring(header) if isinstance(header, str) else header
    row = table.data[0]
    xs, ys = np.asarray(row["X"], float), np.asarray(row["Y"], float)
    values = np.asarray(row["VALUES"], float).reshape(len(ys), len(xs), 2)[..., 0]
    px = hdr["CRVAL1W"] + hdr["CDELT1W"] * (x + 1 - hdr["CRPIX1W"])
    py = hdr["CRVAL2W"] + hdr["CDELT2W"] * (y + 1 - hdr["CRPIX2W"])
    fx = np.interp(px, xs, np.arange(len(xs)))
    fy = np.interp(py, ys, np.arange(len(ys)))
    ix, iy = min(int(fx), len(xs) - 2), min(int(fy), len(ys) - 2)
    tx, ty = fx - ix, fy - iy
    return float((1 - tx) * (1 - ty) * values[iy, ix] + tx * (1 - ty) * values[iy, ix + 1]
                 + (1 - tx) * ty * values[iy + 1, ix] + tx * ty * values[iy + 1, ix + 1])


def wavelength_astropy(header, table, x, y):
    hdr = fits.Header.fromstring(header) if isinstance(header, str) else header
    wcs = WCS(hdr, fobj=fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(header=hdr), table]), key="W")
    wcs.sip = None
    wl, _ = wcs.pixel_to_world(x, y)
    return float(wl.to_value("um"))


# ---------------------------------------------------------------- Gaia DR3


GAIA_ESA = "https://gea.esac.esa.int/tap-server/tap/sync"


def _gaia_esa(where):
    """Gaia DR3 from ESA's own archive (usually a few seconds)."""
    r = _session.post(GAIA_ESA, timeout=30, data={
        "REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv",
        "QUERY": f"SELECT ra, dec, pmra, pmdec, phot_g_mean_mag FROM gaiadr3.gaia_source WHERE {where}"})
    r.raise_for_status()
    rows = list(csv.reader(io.StringIO(r.text)))
    if not rows or rows[0][:2] != ["ra", "dec"]:
        raise IOError("unexpected Gaia answer")
    return [[float(v) if v else 0.0 for v in row] for row in rows[1:]]


def _gaia_irsa(where):
    """Gaia DR3 from IRSA's copy (fallback; slow in crowded fields)."""
    t = _tap(f"SELECT ra, dec, pmra, pmdec, phot_g_mean_mag FROM gaia_dr3_source WHERE {where}")
    cols = [np.ma.filled(np.ma.asarray(t[c], float), 0.0) for c in ("ra", "dec", "pmra", "pmdec", "phot_g_mean_mag")]
    return [list(r) for r in zip(*cols)]


def query_gaia(ra, dec, radius_deg, gmax=19.0):
    """Gaia DR3 stars around (ra, dec): positions at epoch 2016.0, proper motions (mas/yr) and G magnitude."""
    key = hashlib.sha1(f"gaia|{ra:.4f}|{dec:.4f}|{radius_deg:.4f}|{gmax}".encode()).hexdigest()
    path = _cache_path(key + ".npz")
    if os.path.exists(path):
        z = np.load(path)
        return {k: z[k] for k in z.files}
    where = (f"1 = CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, {radius_deg:.5f})) "
             f"AND phot_g_mean_mag < {gmax}")
    try:
        rows = _gaia_esa(where)
    except Exception:
        rows = _gaia_irsa(where)
    a = np.array(rows, float).reshape(-1, 5)
    out = {"ra": a[:, 0], "dec": a[:, 1], "pmra": a[:, 2], "pmdec": a[:, 3], "g": a[:, 4]}
    np.savez_compressed(path, **out)
    return out
