"""End-to-end check on real data: python scripts/smoke_test.py "WISE J085510.83-071442.5" OUT_DIR [fov] [detectors]"""
import os
import sys
import time

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from skyblink import archive, imaging  # noqa: E402

name, out = sys.argv[1], sys.argv[2]
fov = float(sys.argv[3]) if len(sys.argv) > 3 else 4.0
dets = [int(d) for d in sys.argv[4].split(",")] if len(sys.argv) > 4 else [4, 5, 6]
os.makedirs(out, exist_ok=True)

t0 = time.time()
tgt = archive.resolve_name(name)
print("target", tgt)
items = archive.query_images(tgt["ra"], tgt["dec"])
mid_year = imaging.mjd_to_year(np.median([i["mjd"] for i in items]))
ra0, dec0 = imaging.propagate(tgt["ra"], tgt["dec"], tgt["pmra"], tgt["pmdec"], tgt["epoch"], mid_year)
ra0, dec0 = float(ra0), float(dec0)
print(f"{len(items)} images, centre at {mid_year:.2f}: {ra0:.5f} {dec0:.5f}  [{time.time() - t0:.1f}s]")

size = fov / 60 * 1.5
blink_items = imaging.select_for_blink(items, dets)
spec_items = imaging.select_for_spectrum(items)
wanted = list({it["uri"]: it for it in blink_items + spec_items}.values())
cuts = archive.fetch_cutouts(wanted, ra0, dec0, size, progress=lambda d, n: d % 50 == 0 and print(f"  {d}/{n}"))
print(f"{len(cuts)} cutouts  [{time.time() - t0:.1f}s]")

grid, shape = imaging.make_grid(ra0, dec0, fov)
epochs = imaging.build_epochs(cuts, grid, shape, dets)
print("epochs:", [(e["date"], e["n"]) for e in epochs], f"[{time.time() - t0:.1f}s]")
lo, hi = imaging.display_levels([e["image"] for e in epochs])
for k, e in enumerate(epochs):
    Image.fromarray(imaging.to_rgb(e["image"], lo, hi)).resize((shape[1] * 5, shape[0] * 5), Image.LANCZOS).save(
        os.path.join(out, f"epoch{k}.png"))
diffs = imaging.differences(epochs)
if diffs:
    s = imaging.diff_scale(diffs)
    Image.fromarray(imaging.diff_to_rgb(diffs[-1], s)).resize((shape[1] * 5, shape[0] * 5), Image.LANCZOS).save(
        os.path.join(out, "diff_last.png"))
# where the target should be in each epoch
for e in epochs:
    r, d = imaging.propagate(tgt["ra"], tgt["dec"], tgt["pmra"], tgt["pmdec"], tgt["epoch"], e["year"])
    x, y = imaging.to_display_xy(grid, shape, r, d)
    print(f"  {e['date']}: target at display px ({float(x):.1f}, {float(y):.1f})")

waves = {}
for it in items:
    if it["detector"] not in waves:
        waves[it["detector"]] = archive.wave_table(it, ra0, dec0)
spec = imaging.spectrum(cuts, waves, tgt["ra"], tgt["dec"], tgt["pmra"], tgt["pmdec"], tgt["epoch"])
print(f"spectrum: {len(spec['flux'])} points, {spec['masked']} images with a saturated core")
if len(spec["flux"]):
    for lo_w, hi_w in [(0.75, 1.1), (1.1, 1.6), (1.6, 2.4), (2.4, 3.8), (3.8, 4.4), (4.4, 5.0)]:
        m = (spec["wavelength"] >= lo_w) & (spec["wavelength"] < hi_w)
        if m.any():
            print(f"  {lo_w}-{hi_w} um: median {np.median(spec['flux'][m]):.4f} mJy (n={m.sum()})")
print(f"total {time.time() - t0:.1f}s")
