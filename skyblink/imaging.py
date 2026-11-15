"""Turning SPHEREx cutouts into aligned epoch images, difference images and spectra."""
import logging
import re
import warnings

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area
from reproject import reproject_interp

from . import archive

logging.getLogger("astropy").setLevel(logging.ERROR)
warnings.filterwarnings("ignore")

GRID_SCALE = 3.0       # arcsec per display pixel (SPHEREx pixels are ~6.15", so 2x oversampled)
VISIT_GAP = 20.0       # days without data that separate two visits
MAX_VISIT_SPAN = 45.0  # days; longer stretches (deep fields near the ecliptic poles) are split
MAX_PER_EPOCH = 25     # images combined per epoch, chosen evenly in time
MAX_EPOCHS = 8         # deep fields near the ecliptic poles have dozens of visits; keep evenly spaced ones
SPEC_PER_DETECTOR = 45 # images per detector used for a spectrum


def mjd_to_year(mjd):
    return float(Time(mjd, format="mjd").decimalyear)


def mjd_to_date(mjd):
    return Time(mjd, format="mjd").to_datetime().strftime("%Y-%m-%d")


def propagate(ra, dec, pmra, pmdec, epoch, year):
    """Move a position by its proper motion (mas/yr) from `epoch` to `year` (decimal years)."""
    dt = np.asarray(year, float) - epoch
    ra = np.asarray(ra, float) + np.asarray(pmra, float) / 3.6e6 * dt / np.cos(np.radians(dec))
    dec = np.asarray(dec, float) + np.asarray(pmdec, float) / 3.6e6 * dt
    return ra, dec


def make_grid(ra, dec, fov_arcmin):
    """A north-up TAN grid centred on the target."""
    n = int(round(fov_arcmin * 60 / GRID_SCALE)) | 1
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [ra, dec]
    w.wcs.crpix = [(n + 1) / 2, (n + 1) / 2]
    w.wcs.cdelt = [-GRID_SCALE / 3600, GRID_SCALE / 3600]
    return w, (n, n)


def epoch_groups(mjds):
    """Split observation times into survey visits: index lists, oldest first."""
    order = np.argsort(mjds)
    groups, cur = [], [order[0]]
    for i in order[1:]:
        if mjds[i] - mjds[cur[-1]] > VISIT_GAP or mjds[i] - mjds[cur[0]] > MAX_VISIT_SPAN:
            groups.append(cur)
            cur = [i]
        else:
            cur.append(i)
    groups.append(cur)
    if len(groups) > MAX_EPOCHS:
        groups = [groups[int(round(k))] for k in np.linspace(0, len(groups) - 1, MAX_EPOCHS)]
    out = []
    for g in groups:
        if len(g) > MAX_PER_EPOCH:
            g = [g[int(k)] for k in np.linspace(0, len(g) - 1, MAX_PER_EPOCH)]
        out.append(g)
    return out


def select_for_blink(items, detectors):
    """The images worth downloading for the time-lapse: chosen detectors, at most MAX_PER_EPOCH per visit."""
    sel = [it for it in items if it["detector"] in detectors]
    if not sel:
        return []
    return [sel[i] for g in epoch_groups(np.array([it["mjd"] for it in sel])) for i in g]


# alternate-WCS keywords such as CRPIX1W, PS1_0W, RADESYSA (but not RADESYS, NAXIS or SIP's A_1_2)
_ALT_KEY = re.compile(r"^([A-Z]+\d+(_\d+)?|WCSAXES|WCSNAME|RADESYS|EQUINOX|LONPOLE|LATPOLE)[A-Z]$")
_wcs_cache = {}


def cut_wcs(cut):
    """Celestial WCS of a cutout, cached: astropy's WCS() is the slowest step of every loop here."""
    key = hash(cut["header"])
    w = _wcs_cache.get(key)
    if w is None:
        hdr = fits.Header.fromstring(cut["header"])
        for k in [k for k in hdr if _ALT_KEY.match(k)]:
            del hdr[k]
        w = _wcs_cache[key] = WCS(hdr)
        if len(_wcs_cache) > 5000:
            _wcs_cache.clear()
    return w


def select_for_spectrum(items):
    """At most SPEC_PER_DETECTOR images per detector, spread evenly in time (they sample different wavelengths)."""
    out = []
    for det in sorted({it["detector"] for it in items}):
        sel = sorted((it for it in items if it["detector"] == det), key=lambda it: it["mjd"])
        if len(sel) > SPEC_PER_DETECTOR:
            sel = [sel[int(k)] for k in np.linspace(0, len(sel) - 1, SPEC_PER_DETECTOR)]
        out += sel
    return out


def _reproject(cut, grid, shape):
    arr, foot = reproject_interp((cut["data"], cut_wcs(cut)), grid, shape_out=shape, order="bilinear")
    arr[foot < 0.5] = np.nan
    return arr


def fill_holes(img, iters=15):
    """Fill empty pixels with their brightest neighbour, repeatedly.

    The cores of very bright stars are empty in every image (the detector could not measure them);
    this paints them white instead of leaving a black hole in the picture. Display only.
    """
    a = img.copy()
    h, w = a.shape
    holes = ~np.isfinite(a)
    if not holes.any() or holes.all():
        return a
    bright, peak = np.nanpercentile(a, 99), np.nanmax(a)
    for _ in range(iters):
        empty = ~np.isfinite(a)
        if not empty.any():
            break
        p = np.pad(a, 1, constant_values=np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m = np.nanmax(np.stack([p[dy:dy + h, dx:dx + w] for dy in range(3) for dx in range(3)]), axis=0)
        a[empty] = m[empty]
    # holes inside a bright star become one even white disc; holes at the image edge keep the background
    a[holes & (a > bright)] = peak
    return a


def build_epochs(cuts, grid, shape, detectors):
    """Aligned median image per visit, using only the chosen detectors."""
    sel = [c for c in cuts if c["detector"] in detectors]
    if not sel:
        return []
    epochs = []
    for g in epoch_groups(np.array([c["mjd"] for c in sel])):
        frames = []
        for i in g:
            a = _reproject(sel[i], grid, shape)
            if np.isfinite(a).mean() < 0.3:
                continue
            frames.append(a - np.nanmedian(a))
        if not frames:
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            img = fill_holes(np.nanmedian(np.stack(frames), axis=0))
        mjd = float(np.mean([sel[i]["mjd"] for i in g]))
        epochs.append({"image": img, "mjd": mjd, "year": mjd_to_year(mjd), "date": mjd_to_date(mjd), "n": len(frames)})
    return epochs


# ---------------------------------------------------------------- display


LUT_STOPS = [(0.0, (4, 7, 18)), (0.35, (36, 64, 128)), (0.62, (92, 200, 232)), (0.85, (255, 214, 160)), (1.0, (255, 255, 255))]
_x = np.linspace(0, 1, 256)
LUT = np.stack([np.interp(_x, [s for s, _ in LUT_STOPS], [c[i] for _, c in LUT_STOPS]) for i in range(3)], 1).astype(np.uint8)


def display_levels(images):
    vals = np.concatenate([i[np.isfinite(i)] for i in images])
    return float(np.percentile(vals, 20)), float(np.percentile(vals, 99.3))


def to_rgb(img, lo, hi, stretch=12.0):
    """Asinh stretch through the SkyBlink colour map."""
    x = np.clip((np.nan_to_num(img, nan=lo) - lo) / max(hi - lo, 1e-9), 0, 1)
    y = np.arcsinh(x * stretch) / np.arcsinh(stretch)
    return LUT[(y * 255).astype(np.uint8)][::-1]  # FITS row 0 is the bottom


def differences(epochs):
    """Each later visit minus the first one."""
    return [e["image"] - epochs[0]["image"] for e in epochs[1:]]


def diff_to_rgb(d, scale):
    x = np.clip(0.5 + np.nan_to_num(d) / (2 * max(scale, 1e-9)), 0, 1)
    g = (x * 255).astype(np.uint8)
    return np.dstack([g, g, g])[::-1]


def diff_scale(diffs):
    vals = np.concatenate([np.abs(d[np.isfinite(d)]) for d in diffs])
    return float(np.percentile(vals, 99.7))


def to_display_xy(grid, shape, ra, dec):
    """Grid pixel coordinates with y flipped to match the displayed (north-up) image."""
    x, y = grid.world_to_pixel_values(ra, dec)
    return np.asarray(x), shape[0] - 1 - np.asarray(y)


# ---------------------------------------------------------------- sources and spectra


def find_peaks(img, nsigma=6.0, max_n=25):
    """Bright local maxima in an image: (x, y, value), brightest first."""
    a = np.nan_to_num(img, nan=np.nanmedian(img))
    med = np.median(a)
    sig = 1.4826 * np.median(np.abs(a - med))
    core = a[2:-2, 2:-2]
    is_max = np.ones_like(core, bool)
    for dy in (-2, -1, 0, 1, 2):
        for dx in (-2, -1, 0, 1, 2):
            if dx or dy:
                is_max &= core >= a[2 + dy:a.shape[0] - 2 + dy, 2 + dx:a.shape[1] - 2 + dx]
    ys, xs = np.nonzero(is_max & (core > med + nsigma * sig))
    vals = core[ys, xs]
    order = np.argsort(vals)[::-1][:max_n]
    return [(int(xs[i] + 2), int(ys[i] + 2), float(vals[i])) for i in order]


def spectrum(cuts, waves, ra, dec, pmra=0.0, pmdec=0.0, epoch=2000.0):
    """Aperture photometry of one source in every cutout: wavelength (um), flux (mJy), error, detector, mjd.

    `masked` counts images whose aperture was mostly flagged: for bright stars that means a saturated core.
    """
    rows, masked = [], 0
    yy, xx = np.mgrid[-7:8, -7:8]
    rr = np.hypot(xx, yy)
    ann = (rr >= 4) & (rr <= 7)
    for c in cuts:
        table = waves.get(c["detector"])
        if table is None:
            continue
        w = cut_wcs(c)
        r_t, d_t = propagate(ra, dec, pmra, pmdec, epoch, mjd_to_year(c["mjd"]))
        x, y = w.world_to_pixel_values(r_t, d_t)
        d = np.where(c["soft"], np.nan, c["data"])
        xi, yi = int(round(float(x))), int(round(float(y)))
        if not (7 <= xi < d.shape[1] - 7 and 7 <= yi < d.shape[0] - 7):
            continue
        box = d[yi - 7:yi + 8, xi - 7:xi + 8]
        sky = box[ann]
        if np.isfinite(sky).sum() < 20:
            continue
        bkg, noise = np.nanmedian(sky), 1.4826 * np.nanmedian(np.abs(sky - np.nanmedian(sky)))
        # forced photometry at the predicted position; re-centre only on a clear detection,
        # otherwise picking the brightest pixel would turn noise into fake flux
        inner = np.where(rr <= 1.5, np.nan_to_num(box, nan=-np.inf), -np.inf)
        py, px = np.unravel_index(np.argmax(inner), inner.shape)
        cy, cx = (yi - 7 + py, xi - 7 + px) if box[py, px] - bkg > 5 * noise else (yi, xi)
        ap = d[cy - 1:cy + 2, cx - 1:cx + 2]
        if np.isfinite(ap).sum() < 7:
            masked += 1
            continue
        flux = np.nanmean(ap - bkg) * 9
        omega = proj_plane_pixel_area(w) * (np.pi / 180) ** 2  # sr per pixel
        to_mjy = omega * 1e9  # MJy/sr * sr = MJy; 1 MJy = 1e9 mJy
        rows.append((archive.wavelength(c["header"], table, cx, cy), flux * to_mjy, noise * 3 * to_mjy, c["detector"], c["mjd"]))
    a = np.array(rows) if rows else np.zeros((0, 5))
    order = np.argsort(a[:, 0])
    return {"wavelength": a[order, 0], "flux": a[order, 1], "error": a[order, 2],
            "detector": a[order, 3].astype(int), "mjd": a[order, 4], "masked": masked}
