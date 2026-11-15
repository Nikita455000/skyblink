"""Pack the cutouts of each challenge target's default view (4' field, detectors D4-D6, plus the spectrum
images) into skyblink/data/seed/<target>.npz, so the deployed app opens them in seconds.

Mirrors what app.py does for a search, so the cutout keys match. Needs the target's saved image list
in skyblink/data/targets/ (see save_target_lists.py).

Usage: python scripts/build_seed.py [target ...]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from skyblink import archive, imaging  # noqa: E402

NAMES = {"wise_0855": "WISE J085510.83-071442.5", "barnard": "Barnard's Star", "eps_ind_b": "eps Ind B",
         "rw_aur": "RW Aur", "r_dra": "R Dra"}
FOV, DETS = 4.0, [4, 5, 6]
OUT = archive.SEED_DIR
TARGETS = os.path.join(os.path.dirname(archive.SEED_DIR), "targets")

os.makedirs(OUT, exist_ok=True)
for key in sys.argv[1:] or NAMES:
    if not os.path.exists(os.path.join(TARGETS, key + ".json")):
        print(key, "- no saved image list yet, skipped")
        continue
    t = archive.resolve_name(NAMES[key])
    ra_now, dec_now = (float(v) for v in imaging.propagate(t["ra"], t["dec"], t["pmra"], t["pmdec"], t["epoch"], 2026.0))
    items = archive.query_images(round(ra_now, 5), round(dec_now, 5))
    mid = imaging.mjd_to_year(float(np.median([i["mjd"] for i in items])))
    ra0, dec0 = (float(v) for v in imaging.propagate(t["ra"], t["dec"], t["pmra"], t["pmdec"], t["epoch"], mid))
    size = imaging.cutout_size(FOV)
    wanted = list({it["uri"]: it for it in imaging.select_for_blink(items, DETS) + imaging.select_for_spectrum(items)}.values())
    cuts = archive.fetch_cutouts(wanted, ra0, dec0, size, progress=lambda d, n: d % 100 == 0 and print(f"  {d}/{n}"))
    arrays = {}
    by_uri = {c["uri"]: c for c in cuts}
    for it in wanted:
        c = by_uri.get(it["uri"])
        if c is None:
            continue
        k = archive.cutout_key(it, ra0, dec0, size)
        arrays[k + "_data"], arrays[k + "_soft"], arrays[k + "_hdr"] = c["data"], c["soft"], np.array(c["header"])
    path = os.path.join(OUT, key + ".npz")
    np.savez_compressed(path, **arrays)
    print(f"{key}: {len(arrays) // 3} of {len(wanted)} cutouts, {os.path.getsize(path) / 1e6:.1f} MB")
