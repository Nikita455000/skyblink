"""SkyBlink: find what moves in the infrared sky. NASA Space Apps 2026, team DeepSpace.

Run locally with:  streamlit run app.py
"""
import concurrent.futures
import io
import json

import numpy as np
import plotly.graph_objects as go
import streamlit as st
import streamlit.components.v1 as components
from astropy import units as u
from astropy.coordinates import SkyCoord
from PIL import Image, ImageDraw

from skyblink import archive, imaging, viewer

st.set_page_config(page_title="SkyBlink — find what moves", page_icon="🔭", layout="wide")

TARGETS = {
    "Barnard’s Star": "Barnard's Star",
    "Epsilon Indi B": "eps Ind B",
    "WISE 0855": "WISE J085510.83-071442.5",
    "RW Aurigae": "RW Aur",
    "R Draconis": "R Dra",
}
DET_LABELS = {d: f"D{d} · {lo:.2f}–{hi:.2f} µm" for d, (lo, hi) in archive.DETECTORS.items()}
DET_COLORS = {1: "#3B6FD8", 2: "#4FA3E0", 3: "#5CC8E8", 4: "#E8C547", 5: "#FF9A5A", 6: "#FF7A45"}
DEFAULT_DETS = [4, 5, 6]

st.markdown("""<style>
  .block-container { padding-top: 2rem; }
  h1, h2, h3 { letter-spacing: -0.01em; }
  .sb-eyebrow { color:#FF7A45; font-weight:700; letter-spacing:.18em; text-transform:uppercase; font-size:.8rem; }
  .sb-muted { color:#B8C4D6; }
</style>""", unsafe_allow_html=True)


# ---------------------------------------------------------------- cached archive calls


@st.cache_data(show_spinner=False, ttl=86400)
def resolve(name):
    return archive.resolve_name(name)


@st.cache_data(show_spinner=False, ttl=6 * 3600)
def images_at(ra, dec):
    return archive.query_images(ra, dec)


@st.cache_resource
def _gaia_pool():
    return concurrent.futures.ThreadPoolExecutor(max_workers=4)


GAIA_POOL = _gaia_pool()


def parse_coords(text):
    text = text.strip().replace(",", " ")
    unit = (u.hourangle, u.deg) if any(ch in text for ch in "h:") else (u.deg, u.deg)
    c = SkyCoord(text, unit=unit, frame="icrs")
    return float(c.ra.deg), float(c.dec.deg)


# ---------------------------------------------------------------- search inputs

qp = st.query_params
if "req" not in st.session_state and ("target" in qp or "coords" in qp):
    st.session_state.req = {
        "name": qp.get("target", ""),
        "coords": qp.get("coords", ""),
        "fov": float(qp.get("fov", 4.0)),
        "dets": [int(d) for d in qp.get("det", "456") if d in "123456"] or DEFAULT_DETS,
    }

with st.sidebar:
    st.markdown("## 🔭 SkyBlink")
    st.caption("Find what moves in the infrared sky")
    st.markdown("**Challenge targets**")
    cols = st.columns(2)
    for i, label in enumerate(TARGETS):
        if cols[i % 2].button(label, use_container_width=True):
            st.session_state.req = {"name": TARGETS[label], "coords": "",
                                    "fov": st.session_state.get("fov_in", 4.0),
                                    "dets": st.session_state.get("dets_in", DEFAULT_DETS)}
    req = st.session_state.get("req", {})
    with st.form("search"):
        name = st.text_input("Object name", value=req.get("name", ""), placeholder="WISE 0855, Barnard's Star, R Dra…")
        coords = st.text_input("…or coordinates (RA Dec)", value=req.get("coords", ""),
                               placeholder="08h55m10.8s -07d14m42.5s  or  133.79 -7.24")
        fov = st.slider("Field of view (arcmin)", 2.0, 10.0, float(req.get("fov", 4.0)), 0.5, key="fov_in")
        dets = st.multiselect("Detectors for the time-lapse", list(archive.DETECTORS), default=req.get("dets", DEFAULT_DETS),
                              format_func=DET_LABELS.get, key="dets_in")
        if st.form_submit_button("Search the sky", type="primary", use_container_width=True):
            st.session_state.req = {"name": name, "coords": coords, "fov": fov, "dets": sorted(dets) or DEFAULT_DETS}
    st.divider()
    st.caption("Data: NASA SPHEREx and ESA Gaia DR3 via NASA/IPAC IRSA · names via CDS Sesame")
    st.caption("NASA Space Apps 2026 · Planet X and SPHEREx · Team DeepSpace")


# ---------------------------------------------------------------- the pipeline


def run_search(req):
    key = json.dumps(req, sort_keys=True)
    if st.session_state.get("result_key") == key:
        return st.session_state.result
    with st.status("Searching the SPHEREx archive…", expanded=True) as status:
        if req["coords"].strip():
            try:
                ra, dec = parse_coords(req["coords"])
            except Exception:
                status.update(label="Could not read those coordinates", state="error")
                return None
            tgt = {"name": req["coords"].strip(), "ra": ra, "dec": dec, "pmra": 0.0, "pmdec": 0.0, "epoch": 2000.0}
        else:
            status.write(f"Looking up **{req['name']}**…")
            try:
                tgt = resolve(req["name"]) if req["name"].strip() else None
            except Exception:
                status.update(label="The name resolver is not answering — try coordinates instead", state="error")
                return None
            if tgt is None:
                status.update(label=f"No object called “{req['name']}” was found", state="error")
                return None
        # high proper-motion stars move arcminutes in decades: search where the object is now
        ra_now, dec_now = (float(v) for v in imaging.propagate(tgt["ra"], tgt["dec"], tgt["pmra"], tgt["pmdec"], tgt["epoch"], 2026.0))
        status.write("Finding SPHEREx images of this spot…")
        try:
            items = images_at(round(ra_now, 5), round(dec_now, 5))
        except archive.ArchiveUnavailable:
            status.update(label="NASA’s IRSA search service is not answering right now", state="error")
            st.warning("IRSA (NASA/IPAC Infrared Science Archive) is temporarily not answering image searches. "
                       "The five challenge targets on the left still work from saved image lists — "
                       "or try this position again in a few minutes.")
            return None
        if not items:
            status.update(label="SPHEREx has not observed this position yet", state="error")
            return None
        mid = imaging.mjd_to_year(float(np.median([i["mjd"] for i in items])))
        ra0, dec0 = (float(v) for v in imaging.propagate(tgt["ra"], tgt["dec"], tgt["pmra"], tgt["pmdec"], tgt["epoch"], mid))
        sel = imaging.select_for_blink(items, req["dets"])
        size = req["fov"] / 60 * 1.5
        # Gaia is optional: ask for it in the background while the images download
        gaia_job = GAIA_POOL.submit(archive.query_gaia, round(ra0, 4), round(dec0, 4), round(req["fov"] / 60 * 0.75, 4))
        status.write(f"Downloading {len(sel)} image cutouts from IRSA (cached after the first time)…")
        bar = st.progress(0.0)
        cuts = archive.fetch_cutouts(sel, ra0, dec0, size, progress=lambda d, n: bar.progress(d / n, text=f"{d} / {n} images"))
        status.write("Aligning the images and building the time-lapse…")
        grid, shape = imaging.make_grid(ra0, dec0, req["fov"])
        epochs = imaging.build_epochs(cuts, grid, shape, req["dets"])
        if not epochs:
            status.update(label="No usable images for these detectors here", state="error")
            return None
        status.write("Adding Gaia DR3 stars…")
        try:
            gaia = gaia_job.result(timeout=20)
        except Exception:
            gaia = None
            status.write("Gaia DR3 is not answering right now — showing the images without Gaia circles.")
        status.update(label=f"Ready — {len(epochs)} visits of {tgt['name']}", state="complete", expanded=False)
    res = {"req": req, "tgt": tgt, "items": items, "ra0": ra0, "dec0": dec0, "size": size, "mid": mid,
           "grid": grid, "shape": shape, "epochs": epochs, "gaia": gaia}
    st.session_state.result_key, st.session_state.result = key, res
    st.session_state.pop("spectrum", None)
    return res


def target_xy(res, year):
    t = res["tgt"]
    r, d = imaging.propagate(t["ra"], t["dec"], t["pmra"], t["pmdec"], t["epoch"], year)
    x, y = imaging.to_display_xy(res["grid"], res["shape"], r, d)
    return [round(float(x), 2), round(float(y), 2)]


def gaia_xy(res, year):
    g = res["gaia"]
    if g is None or not len(g["ra"]):
        return np.zeros(0), np.zeros(0), np.zeros(0, bool)
    r, d = imaging.propagate(g["ra"], g["dec"], g["pmra"], g["pmdec"], 2016.0, year)
    x, y = imaging.to_display_xy(res["grid"], res["shape"], r, d)
    n = res["shape"][0]
    return x, y, (x > -1) & (x < n) & (y > -1) & (y < n)


def gaia_list(res, year, max_n=40):
    """Gaia stars inside the frame, brightest first, capped so crowded fields stay readable."""
    x, y, inside = gaia_xy(res, year)
    idx = np.nonzero(inside)[0]
    idx = idx[np.argsort(res["gaia"]["g"][idx])][:max_n] if len(idx) else idx
    return [[round(float(x[i]), 2), round(float(y[i]), 2)] for i in idx]


def viewer_payload(res):
    ep = res["epochs"]
    lo, hi = imaging.display_levels([e["image"] for e in ep])
    blink = [{"src": viewer.png_data_url(imaging.to_rgb(e["image"], lo, hi)),
              "label": f"Visit {i + 1} of {len(ep)} · {e['date']} · {e['n']} images",
              "target": target_xy(res, e["year"]), "gaia": gaia_list(res, e["year"])} for i, e in enumerate(ep)]
    diffs = imaging.differences(ep)
    diff = []
    if diffs:
        s = imaging.diff_scale(diffs)
        diff = [{"src": viewer.png_data_url(imaging.diff_to_rgb(d, s)),
                 "label": f"Visit {i + 2} minus visit 1 · white = now, black = then",
                 "target": target_xy(res, ep[i + 1]["year"]), "gaia": gaia_list(res, ep[i + 1]["year"])}
                for i, d in enumerate(diffs)]
    stack = np.nanmean(np.stack([e["image"] for e in ep]), axis=0)
    stk = [{"src": viewer.png_data_url(imaging.to_rgb(stack, lo, hi)), "label": f"All {len(ep)} visits combined",
            "target": target_xy(res, res["mid"]), "gaia": gaia_list(res, res["mid"])}]
    return {"n": res["shape"][0], "scale": imaging.GRID_SCALE, "maxSize": 680, "blink": blink, "diff": diff, "stack": stk}, stack, (lo, hi)


def gif_bytes(res, lo, hi):
    frames = []
    for e in res["epochs"]:
        im = Image.fromarray(imaging.to_rgb(e["image"], lo, hi)).resize((480, 480), Image.BICUBIC)
        x, y = target_xy(res, e["year"])
        k = 480 / res["shape"][0]
        cx, cy = (x + 0.5) * k, (y + 0.5) * k
        d = ImageDraw.Draw(im)
        d.ellipse((cx - 18, cy - 18, cx + 18, cy + 18), outline=(255, 122, 69), width=3)
        d.text((12, 10), f"SkyBlink · {res['tgt']['name']} · {e['date']}", fill=(243, 241, 234))
        frames.append(im.convert("P", palette=Image.ADAPTIVE))
    buf = io.BytesIO()
    frames[0].save(buf, format="GIF", save_all=True, append_images=frames[1:], duration=600, loop=0)
    return buf.getvalue()


def share_params(req):
    p = {"fov": f"{req['fov']:g}", "det": "".join(str(d) for d in req["dets"])}
    if req["coords"].strip():
        p["coords"] = req["coords"].strip()
    else:
        p["target"] = req["name"]
    return p


# ---------------------------------------------------------------- page

req = st.session_state.get("req")
if not req:
    st.markdown('<div class="sb-eyebrow">NASA SPHEREx · Planet X challenge</div>', unsafe_allow_html=True)
    st.title("Find what moves in the infrared sky")
    st.markdown(
        "SPHEREx maps the **whole sky every six months in 102 infrared colors**. SkyBlink downloads every "
        "SPHEREx image of the spot you choose, lines them up and plays them as a time-lapse — so comets, "
        "asteroids, brown dwarfs and variable stars jump out.")
    c1, c2, c3 = st.columns(3)
    c1.markdown("#### 1 · Search\nType a name or coordinates, or pick a challenge target on the left.")
    c2.markdown("#### 2 · Watch\n**Blink** plays the visits; **Difference** hides everything that stayed still.")
    c3.markdown("#### 3 · Check\nClick any source to measure its spectrum; cyan circles mark known Gaia stars.")
    st.info("👈 Start with **WISE 0855** — the coldest known brown dwarf, crawling across the sky.")
    st.stop()

res = run_search(req)
if res is None:
    st.stop()
st.query_params.from_dict(share_params(req))

tgt, ep = res["tgt"], res["epochs"]
st.markdown('<div class="sb-eyebrow">SkyBlink · SPHEREx time-lapse</div>', unsafe_allow_html=True)
st.title(tgt["name"])
pos = SkyCoord(res["ra0"] * u.deg, res["dec0"] * u.deg).to_string("hmsdms", precision=1)
pm = np.hypot(tgt["pmra"], tgt["pmdec"]) / 1000
st.markdown(f'<span class="sb-muted">Position in {res["mid"]:.1f}: {pos} · proper motion {pm:.2f}″/yr · '
            f'detectors {", ".join(f"D{d}" for d in req["dets"])}</span>', unsafe_allow_html=True)

payload, stack, (lo, hi) = viewer_payload(res)
left, right = st.columns([3, 2], gap="large")
with left:
    components.html(viewer.viewer_html(payload), height=800)
with right:
    span = (ep[-1]["year"] - ep[0]["year"]) if len(ep) > 1 else 0
    m1, m2, m3 = st.columns(3)
    m1.metric("Visits", len(ep))
    m2.metric("Images", sum(e["n"] for e in ep))
    m3.metric("Time span", f"{span:.1f} yr")
    expected = pm * span / 6.15
    if len(ep) < 2:
        st.warning("SPHEREx has visited this spot only once so far, so there is no motion to see yet.")
    elif pm > 0.3:
        st.success(f"Expected shift over this span: about **{pm * span:.1f}″** "
                   f"(≈ {expected:.1f} SPHEREx pixels). Watch the orange crosshair in **Difference** mode.")
    st.markdown("**How to read it**")
    st.markdown("- **Blink** — each frame is one SPHEREx visit (median of its images).\n"
                "- **Difference** — later visit minus the first: a white–black pair means something moved; "
                "a lone white or black dot means it got brighter or fainter.\n"
                "- **Cyan circles** — known Gaia DR3 stars at their predicted position. "
                "Something moving with no circle is worth a closer look.")
    st.dataframe({"Visit": [i + 1 for i in range(len(ep))], "Date": [e["date"] for e in ep],
                  "Images": [e["n"] for e in ep]}, hide_index=True, use_container_width=True)
    st.download_button("Download time-lapse (GIF)", gif_bytes(res, lo, hi),
                       file_name=f"skyblink_{tgt['name'].replace(' ', '_')}.gif", mime="image/gif",
                       use_container_width=True)
    st.caption("🔗 The page address now holds this search — copy it to share this exact view.")

# ---------------------------------------------------------------- spectrum

st.divider()
st.subheader("Spectrum · 102 colors")
st.markdown('<span class="sb-muted">Click a source on the map (orange = target, cyan = Gaia star, white = bright source) '
            'and measure its spectrum across 0.75–5 µm.</span>', unsafe_allow_html=True)

n = res["shape"][0]
pick_fig = go.Figure(go.Image(z=imaging.to_rgb(stack, lo, hi)))
tx, ty = target_xy(res, res["mid"])
gx, gy, inside = gaia_xy(res, res["mid"])
g = res["gaia"]
if g is not None and inside.any():
    idx = np.nonzero(inside)[0]
    idx = idx[np.argsort(g["g"][idx])][:40]
    pick_fig.add_trace(go.Scatter(
        x=gx[idx], y=gy[idx], mode="markers", name="Gaia DR3",
        marker=dict(size=14, color="rgba(0,0,0,0)", line=dict(color="#5CC8E8", width=2)),
        customdata=[["gaia", g["ra"][i], g["dec"][i], g["pmra"][i], g["pmdec"][i], f"Gaia star G={g['g'][i]:.1f}"] for i in idx],
        hovertemplate="%{customdata[5]}<extra></extra>"))
peaks = imaging.find_peaks(stack)
if peaks:
    pr, pd = res["grid"].pixel_to_world_values([p[0] for p in peaks], [n - 1 - p[1] for p in peaks])
    pick_fig.add_trace(go.Scatter(
        x=[p[0] for p in peaks], y=[p[1] for p in peaks], mode="markers", name="Bright sources",
        marker=dict(size=10, color="rgba(0,0,0,0)", line=dict(color="#F3F1EA", width=1.5)),
        customdata=[["peak", float(a), float(b), 0.0, 0.0, f"Source at {a:.4f}, {b:+.4f}"] for a, b in zip(pr, pd)],
        hovertemplate="%{customdata[5]}<extra></extra>"))
pick_fig.add_trace(go.Scatter(
    x=[tx], y=[ty], mode="markers", name="Target",
    marker=dict(size=20, symbol="circle-open", color="#FF7A45", line=dict(width=3)),
    customdata=[["target", tgt["ra"], tgt["dec"], tgt["pmra"], tgt["pmdec"], tgt["name"]]],
    hovertemplate="%{customdata[5]}<extra></extra>"))
pick_fig.update_layout(height=520, margin=dict(l=0, r=0, t=0, b=0), paper_bgcolor="rgba(0,0,0,0)",
                       plot_bgcolor="rgba(0,0,0,0)", showlegend=True, clickmode="event+select",
                       legend=dict(orientation="h", y=-0.02, font=dict(color="#B8C4D6")))
pick_fig.update_xaxes(visible=False, range=[-0.5, n - 0.5])
pick_fig.update_yaxes(visible=False, range=[n - 0.5, -0.5], scaleanchor="x")

c_map, c_spec = st.columns([2, 3], gap="large")
with c_map:
    event = st.plotly_chart(pick_fig, on_select="rerun", selection_mode="points", key="picker", use_container_width=True)
points = event.selection.points if event and event.selection else []
src = points[0]["customdata"] if points and points[0].get("customdata") else ["target", tgt["ra"], tgt["dec"], tgt["pmra"], tgt["pmdec"], tgt["name"]]
kind, s_ra, s_dec, s_pmra, s_pmdec, s_name = src
s_epoch = {"target": tgt["epoch"], "gaia": 2016.0}.get(kind, res["mid"])

with c_spec:
    st.markdown(f"**Selected:** {s_name}")
    spec_key = f"{s_ra:.5f}|{s_dec:.5f}|{res['ra0']:.4f}"
    spectra = st.session_state.setdefault("spectrum", {})
    spec_items = imaging.select_for_spectrum(res["items"])
    if spec_key not in spectra:
        st.caption(f"Uses {len(spec_items)} SPHEREx images of this spot, spread over all six detectors — "
                   "the first run takes a minute or two.")
        if st.button("Measure spectrum", type="primary"):
            bar = st.progress(0.0)
            cuts = archive.fetch_cutouts(spec_items, res["ra0"], res["dec0"], res["size"],
                                         progress=lambda d, n_: bar.progress(d / n_, text=f"{d} / {n_} images"))
            waves = {}
            for it in spec_items:
                if it["detector"] not in waves:
                    waves[it["detector"]] = archive.wave_table(it, res["ra0"], res["dec0"])
            spectra[spec_key] = imaging.spectrum(cuts, waves, s_ra, s_dec, s_pmra, s_pmdec, s_epoch)
            st.rerun()
    else:
        sp = spectra[spec_key]
        if not len(sp["flux"]):
            if sp["masked"]:
                st.warning(f"Too bright for SPHEREx: the centre of this source is saturated and flagged in "
                           f"{sp['masked']} images, so its spectrum can’t be measured. Try a fainter source.")
            else:
                st.warning("This source falls too close to the edge of the images to measure.")
        else:
            fig = go.Figure()
            for d in sorted(set(sp["detector"])):
                m = sp["detector"] == d
                fig.add_trace(go.Scatter(
                    x=sp["wavelength"][m], y=sp["flux"][m], mode="markers", name=f"D{d}",
                    marker=dict(size=7, color=DET_COLORS[d]),
                    error_y=dict(type="data", array=sp["error"][m], thickness=1, width=0, color=DET_COLORS[d])))
            fig.update_layout(height=470, margin=dict(l=10, r=10, t=10, b=10), paper_bgcolor="rgba(0,0,0,0)",
                              plot_bgcolor="#0E1B33", font=dict(color="#F3F1EA"),
                              legend=dict(orientation="h", y=1.08),
                              xaxis=dict(title="Wavelength (µm)", gridcolor="#24365A"),
                              yaxis=dict(title="Flux density (mJy)", gridcolor="#24365A"))
            st.plotly_chart(fig, use_container_width=True)
            skipped = f" {sp['masked']} images skipped because the core was saturated." if sp["masked"] else ""
            st.caption(f"{len(sp['flux'])} measurements: 3×3-pixel aperture photometry on each SPHEREx image, "
                       f"local background removed. SPHEREx pixels are 6″ wide, so close neighbours can blend.{skipped}")

with st.expander("About SkyBlink and the data"):
    st.markdown(
        "**How it works.** SkyBlink asks NASA/IPAC's Infrared Science Archive (IRSA) for every SPHEREx spectral "
        "image that covers your position, downloads small cutouts, removes zodiacal light and flagged pixels, "
        "aligns everything on one north-up grid and combines the images of each visit. Difference images subtract "
        "the first visit from later ones. Positions of high proper-motion stars are propagated from their catalog "
        "epoch, and Gaia DR3 stars are moved to the date of each visit.\n\n"
        "**Data.** NASA SPHEREx spectral images (IRSA) · ESA Gaia DR3 (IRSA copy) · object names from CDS Sesame.\n\n"
        "**Team DeepSpace** — Nikita Tkachenko (team lead), Dmytro Simonov, Yehor Kaminskyi. "
        "NASA Space Apps Challenge 2026, challenge “Planet X and SPHEREx”.")
