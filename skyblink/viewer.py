"""The blink viewer: a small self-contained HTML/JS page embedded in Streamlit."""
import base64
import io
import json

from PIL import Image


def png_data_url(rgb):
    buf = io.BytesIO()
    Image.fromarray(rgb).save(buf, "PNG")
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


TEMPLATE = r"""<!doctype html>
<html><head><meta charset="utf-8">
<style>
  :root { --bg:#0B1426; --panel:#13213A; --line:#34486E; --text:#F3F1EA; --muted:#B8C4D6; --orange:#FF7A45; --cyan:#5CC8E8; }
  html, body { margin:0; background:transparent; color:var(--text); font:14px/1.4 system-ui, -apple-system, "Segoe UI", sans-serif; }
  .bar { display:flex; flex-wrap:wrap; gap:10px; align-items:center; margin:0 0 10px; }
  .seg { display:inline-flex; background:var(--panel); border:1px solid var(--line); border-radius:10px; padding:3px; }
  .seg button { background:transparent; color:var(--muted); border:0; padding:7px 14px; border-radius:7px; font:600 14px system-ui, sans-serif; cursor:pointer; }
  .seg button.on { background:var(--cyan); color:var(--bg); }
  .btn { background:var(--orange); color:var(--bg); border:0; border-radius:9px; padding:8px 16px; font:700 14px system-ui, sans-serif; cursor:pointer; min-width:84px; }
  label { color:var(--muted); display:inline-flex; gap:6px; align-items:center; cursor:pointer; }
  input[type=range] { accent-color:var(--orange); width:110px; }
  input[type=checkbox] { accent-color:var(--cyan); }
  canvas { display:block; border-radius:12px; border:1px solid var(--line); background:#04070f; max-width:100%; }
  .info { display:flex; justify-content:space-between; align-items:center; gap:12px; margin-top:10px; flex-wrap:wrap; }
  #label { font-weight:700; font-size:16px; }
  .dots { display:flex; gap:8px; }
  .dots span { width:12px; height:12px; border-radius:50%; background:var(--line); cursor:pointer; }
  .dots span.on { background:var(--orange); }
</style></head>
<body>
<div class="bar">
  <div class="seg" id="modes">
    <button data-m="blink" class="on">Blink</button><button data-m="diff">Difference</button><button data-m="stack">Stack</button>
  </div>
  <button class="btn" id="play">Pause</button>
  <label>Speed <input type="range" id="speed" min="1" max="8" value="3"></label>
  <label><input type="checkbox" id="gaia" checked> Gaia DR3</label>
  <label><input type="checkbox" id="tgt" checked> Target</label>
  <label><input type="checkbox" id="smooth" checked> Smooth</label>
</div>
<canvas id="cv" width="640" height="640"></canvas>
<div class="info"><span id="label"></span><div class="dots" id="dots"></div></div>
<script>
const D = __DATA__;
const $ = id => document.getElementById(id);
const cv = $('cv'), ctx = cv.getContext('2d');
let mode = 'blink', k = 0, playing = true, timer = null;
const imgs = {};
for (const m of ['blink', 'diff', 'stack']) {
  imgs[m] = D[m].map(f => { const im = new Image(); im.onload = draw; im.src = f.src; return im; });
}
function fit() {
  const w = Math.min(document.body.clientWidth - 2, D.maxSize);
  cv.width = w; cv.height = w; draw();
}
function draw() {
  const list = imgs[mode];
  if (!list.length) return;
  const i = k % list.length, im = list[i], f = D[mode][i], s = cv.width / D.n;
  ctx.imageSmoothingEnabled = $('smooth').checked; ctx.imageSmoothingQuality = 'high';
  ctx.fillStyle = '#04070f'; ctx.fillRect(0, 0, cv.width, cv.height);
  if (im.complete) ctx.drawImage(im, 0, 0, cv.width, cv.height);
  if ($('gaia').checked && f.gaia) {
    ctx.strokeStyle = 'rgba(92,200,232,0.85)'; ctx.lineWidth = 1.5;
    for (const [x, y] of f.gaia) { ctx.beginPath(); ctx.arc((x + 0.5) * s, (y + 0.5) * s, Math.max(7, 2.4 * s), 0, 2 * Math.PI); ctx.stroke(); }
  }
  if ($('tgt').checked && f.target) {
    const cx = (f.target[0] + 0.5) * s, cy = (f.target[1] + 0.5) * s, r = Math.max(12, 4 * s), g = r * 0.45;
    ctx.strokeStyle = '#FF7A45'; ctx.lineWidth = 2.5; ctx.beginPath();
    ctx.moveTo(cx - r - 10, cy); ctx.lineTo(cx - g - 6, cy); ctx.moveTo(cx + g + 6, cy); ctx.lineTo(cx + r + 10, cy);
    ctx.moveTo(cx, cy - r - 10); ctx.lineTo(cx, cy - g - 6); ctx.moveTo(cx, cy + g + 6); ctx.lineTo(cx, cy + r + 10);
    ctx.stroke();
  }
  const bar = 60 / D.scale * s;
  ctx.fillStyle = '#F3F1EA'; ctx.fillRect(16, cv.height - 22, bar, 3);
  ctx.font = '600 14px system-ui, sans-serif'; ctx.fillText('1′', 22 + bar, cv.height - 16);
  ctx.fillText('N ↑   E ←', cv.width - 82, 24);
  $('label').textContent = f.label;
  $('dots').innerHTML = list.length > 1 ? list.map((_, j) => `<span class="${j === i ? 'on' : ''}" data-j="${j}"></span>`).join('') : '';
}
function restart() {
  clearInterval(timer);
  $('play').textContent = playing ? 'Pause' : 'Play';
  if (playing && imgs[mode].length > 1) timer = setInterval(() => { k = (k + 1) % imgs[mode].length; draw(); }, 1000 / +$('speed').value);
}
$('modes').onclick = e => {
  const m = e.target.dataset.m; if (!m || !D[m].length) return;
  mode = m; k = 0;
  for (const b of $('modes').children) b.classList.toggle('on', b.dataset.m === m);
  draw(); restart();
};
$('play').onclick = () => { playing = !playing; restart(); };
$('speed').oninput = restart;
for (const id of ['gaia', 'tgt', 'smooth']) $(id).onchange = draw;
$('dots').onclick = e => { if (e.target.dataset.j) { k = +e.target.dataset.j; playing = false; restart(); draw(); } };
window.addEventListener('resize', fit);
if (!D.diff.length) $('modes').children[1].disabled = true;
fit(); restart();
</script></body></html>"""


def viewer_html(payload):
    return TEMPLATE.replace("__DATA__", json.dumps(payload))
