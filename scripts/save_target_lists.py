"""Save SPHEREx image lists for the challenge targets into skyblink/data/targets/, so the app keeps
working for them when IRSA's search services are down. Retries until IRSA answers.

Usage: python scripts/save_target_lists.py [minutes_between_tries] [max_hours]
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from skyblink import archive, imaging  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "skyblink", "data", "targets")
NAMES = {"wise_0855": "WISE J085510.83-071442.5", "barnard": "Barnard's Star", "eps_ind_b": "eps Ind B",
         "rw_aur": "RW Aur", "r_dra": "R Dra"}

wait = float(sys.argv[1]) if len(sys.argv) > 1 else 3
deadline = time.time() + 3600 * (float(sys.argv[2]) if len(sys.argv) > 2 else 3)
todo = [k for k in NAMES if not os.path.exists(os.path.join(OUT, k + ".json"))]
while todo and time.time() < deadline:
    for key in list(todo):
        t = archive.resolve_name(NAMES[key])
        # the app searches at the position the object has in 2026
        ra, dec = (round(float(v), 5) for v in imaging.propagate(t["ra"], t["dec"], t["pmra"], t["pmdec"], t["epoch"], 2026.0))
        try:
            cols, rows = archive._tap_rows(f"""
                SELECT a.uri, p.time_bounds_lower, p.energy_bandpassname
                FROM spherex.artifact a JOIN spherex.plane p ON a.planeid = p.planeid
                WHERE 1 = CONTAINS(POINT('ICRS', {ra:.6f}, {dec:.6f}), p.poly)
                ORDER BY p.time_bounds_lower""", tries=1)
        except archive.ArchiveUnavailable as e:
            print(time.strftime("%H:%M"), key, "-", e, flush=True)
            break
        items = [{"uri": r[0], "mjd": float(r[1]), "detector": int(r[2][-1])} for r in rows]
        with open(os.path.join(OUT, key + ".json"), "w", encoding="utf-8") as f:
            json.dump({"name": key, "ra": ra, "dec": dec, "items": items}, f, separators=(",", ":"))
        print(time.strftime("%H:%M"), key, "saved", len(items), "images", flush=True)
        todo.remove(key)
    if todo:
        time.sleep(wait * 60)
print("done" if not todo else f"gave up, still missing: {todo}")
