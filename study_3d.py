from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import html
from pathlib import Path
import re
import socket
import ssl
import tempfile
import threading
import urllib.parse
import urllib.request
import webbrowser

PLUGIN = {
    "name": "show_3d_study_model",
    "description": (
        "Opens an interactive 3D visualization with webcam hand-gesture controls in the browser "
        "for chemistry molecules (e.g., water, benzene, methane, caffeine, glucose) or physics "
        "simulations (e.g., wave_interference, magnetic_helix, electric_field). "
        "Hand controls: pinch and move one hand to rotate, move two hands apart/together to zoom, "
        "and show a peace sign to reset the view. Camera permission is required."
    ),
    "parameters": {
        "type": "OBJECT",
        "properties": {
            "subject": {
                "type": "STRING",
                "enum": ["chemistry", "physics"],
                "description": "Whether the 3D model is for chemistry or physics."
            },
            "topic": {
                "type": "STRING",
                "description": "Name of the molecule (for chemistry) or physics simulation."
            }
        },
        "required": ["subject", "topic"]
    }
}

_PORT = None
_DIR = Path(tempfile.gettempdir()) / "jarvis_3d_lab"
_CACHE_DIR = _DIR / "molecules_cache"
_TRACKER_MODEL = _DIR / "hand_landmarker.task"

_BUILTIN_SDF = {
    "water": """Water
  Offline3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.1173 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.7572   -0.4692 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.7572   -0.4692 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END
$$$$""",
    "methane": """Methane
  Offline3D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6291    0.6291    0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6291   -0.6291    0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6291    0.6291   -0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6291   -0.6291   -0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END
$$$$""",
    "ammonia": """Ammonia
  Offline3D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.1164 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.9397   -0.2716 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.8138   -0.4699   -0.2716 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8138   -0.4699   -0.2716 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
M  END
$$$$""",
    "carbon dioxide": """CarbonDioxide
  Offline3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    1.1600 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000   -1.1600 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  1  3  2  0  0  0  0
M  END
$$$$""",
    "benzene": """Benzene
  Offline3D

 12 12  0  0  0  0  0  0  0  0999 V2000
    1.3950    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6975    1.2081    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6975    1.2081    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3950    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6975   -1.2081    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6975   -1.2081    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4790    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2395    2.1469    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2395    2.1469    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4790    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2395   -2.1469    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2395   -2.1469    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  1  7  1  0  0  0  0
  2  8  1  0  0  0  0
  3  9  1  0  0  0  0
  4 10  1  0  0  0  0
  5 11  1  0  0  0  0
  6 12  1  0  0  0  0
M  END
$$$$"""
}


def _prepare_tracker_model(timeout: int = 12) -> bool:
    if _TRACKER_MODEL.exists() and _TRACKER_MODEL.stat().st_size > 500_000:
        return True
    url = (
        "https://storage.googleapis.com/mediapipe-models/hand_landmarker/"
        "hand_landmarker/float16/1/hand_landmarker.task"
    )
    try:
        _DIR.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=timeout) as response:
            data = response.read()
        if len(data) > 500_000:
            temp = _TRACKER_MODEL.with_suffix(".task.download")
            temp.write_bytes(data)
            temp.replace(_TRACKER_MODEL)
            return True
    except Exception:
        pass
    return False


def _start_server() -> int:
    global _PORT
    if _PORT is not None:
        return _PORT
    _DIR.mkdir(parents=True, exist_ok=True)
    _CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        _, port = s.getsockname()
    handler = partial(SimpleHTTPRequestHandler, directory=str(_DIR))
    srv = ThreadingHTTPServer(("127.0.0.1", port), handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    _PORT = port
    return port


def _get_molecule_sdf(mol: str) -> str:
    key = mol.strip().lower()
    if key in _BUILTIN_SDF:
        return _BUILTIN_SDF.get(key)

    safe_name = re.sub(r"[^a-z0-9_-]+", "_", key)
    cached_file = _CACHE_DIR / f"{safe_name}.sdf"
    if cached_file.exists():
        return cached_file.read_text(encoding="utf-8")

    encoded = urllib.parse.quote(key)
    endpoints = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/SDF?record_type=3d",
        f"https://cactus.nci.nih.gov/chemical/structure/{encoded}/sdf",
    )
    ssl_ctx = ssl.create_default_context()
    last_err = None

    for url in endpoints:
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=12, context=ssl_ctx) as r:
                data = r.read().decode("utf-8", errors="ignore")
                if "V2000" in data or "V3000" in data:
                    cached_file.write_text(data, encoding="utf-8")
                    return data
        except Exception as err:
            last_err = err

    raise RuntimeError(f"Could not fetch '{mol}' from PubChem or NCI ({last_err}). Try water, benzene, methane, ammonia, or carbon dioxide.")


def _hud_script(port: int | None = None) -> str:
    local_model_url = f"http://127.0.0.1:{port}/hand_landmarker.task" if port else ""
    return """
<div id="hud" style="position:fixed;right:16px;bottom:16px;width:250px;background:#0f172a;border:1px solid #06b6d4;border-radius:12px;padding:10px;z-index:9999;color:#e2e8f0;font-family:sans-serif">
  <div style="position:relative;width:100%;height:140px;background:#020617;border-radius:8px;overflow:hidden">
    <video id="cam" autoplay playsinline muted style="position:absolute;inset:0;width:100%;height:100%;object-fit:cover;transform:scaleX(-1)"></video>
    <canvas id="ovl" width="250" height="140" style="position:absolute;inset:0;width:100%;height:100%;transform:scaleX(-1)"></canvas>
  </div>
  <div id="st" style="margin-top:6px;font-size:12px;font-weight:700;color:#22d3ee;text-align:center">LOADING HAND TRACKER...</div>
  <div style="margin-top:4px;font-size:11px;color:#94a3b8">Pinch+Move: Rotate | 2 Hands: Zoom | Peace: Reset</div>
</div>
<script type="module">
const vid = document.getElementById("cam"), cvs = document.getElementById("ovl"), ctx = cvs.getContext("2d"), st = document.getElementById("st");
let lmkr = null, prev = null, prevD = null, lastR = 0, lastT = -1, stream = null;
let smoothPinch = null, smoothPalm = null, smoothZoomDistance = null, lastDetect = 0;
const TRACKER_URL = "https://cdn.jsdelivr.net/npm/@mediapipe/tasks-vision@0.10.21/vision_bundle.mjs";
const WASM_URL = "https://cdn.jsdelivr.net/npm/@mediapipe/tasks-vision@0.10.21/wasm";
const LOCAL_MODEL_URL = """ + repr(local_model_url) + """;
const REMOTE_MODEL_URL = "https://storage.googleapis.com/mediapipe-models/hand_landmarker/hand_landmarker/float16/1/hand_landmarker.task";
const dist = (a, b) => Math.hypot(a.x - b.x, a.y - b.y);
const palmCenter = h => {
  const ids = [0, 5, 9, 13, 17];
  return {
    x: ids.reduce((sum, id) => sum + h[id].x, 0) / ids.length,
    y: ids.reduce((sum, id) => sum + h[id].y, 0) / ids.length
  };
};
const fingerExtended = (h, mcp, pip, tip) =>
  h[tip].y < h[pip].y - 0.012 && h[pip].y < h[mcp].y + 0.02;
async function makeLandmarker() {
  const { FilesetResolver, HandLandmarker } = await import(TRACKER_URL);
  const vis = await FilesetResolver.forVisionTasks(WASM_URL);
  let modelAssetBuffer;
  for (const modelUrl of [LOCAL_MODEL_URL, REMOTE_MODEL_URL].filter(Boolean)) {
    try {
      const response = await fetch(modelUrl, { cache: "force-cache" });
      if (!response.ok) throw new Error("Model request failed: HTTP " + response.status);
      const candidate = await response.arrayBuffer();
      if (candidate.byteLength < 500000) throw new Error("Hand model download was incomplete");
      modelAssetBuffer = candidate;
      break;
    } catch (e) {
      console.warn("Hand model URL unavailable:", modelUrl, e);
    }
  }
  const baseOptions = modelAssetBuffer
    ? { modelAssetBuffer: new Uint8Array(modelAssetBuffer), delegate: "CPU" }
    : { modelAssetPath: REMOTE_MODEL_URL, delegate: "CPU" };
  return HandLandmarker.createFromOptions(vis, {
    baseOptions: {
      ...baseOptions
    },
    runningMode: "VIDEO",
    numHands: 2,
    minHandDetectionConfidence: 0.55,
    minHandPresenceConfidence: 0.55,
    minTrackingConfidence: 0.55
  });
}
async function boot() {
  try {
    if (!navigator.mediaDevices || !navigator.mediaDevices.getUserMedia) {
      throw new Error("Camera API unavailable");
    }
    st.textContent = "STARTING CAMERA...";
    const cameraPromise = navigator.mediaDevices.getUserMedia({
      video: { facingMode: "user", width: { ideal: 480 }, height: { ideal: 270 }, frameRate: { ideal: 24, max: 30 } },
      audio: false
    });
    const trackerPromise = makeLandmarker();

    cameraPromise.then(async result => {
      stream = result;
      vid.srcObject = stream;
      try {
        await vid.play();
      } catch (_) {}
      st.textContent = lmkr ? "HAND CONTROL READY" : "CAMERA READY - LOADING TRACKER...";
      loop();
    }).catch(e => {
      st.textContent = "CAMERA BUSY / BLOCKED";
      st.title = e && e.message ? e.message : "Camera unavailable";
      st.style.color = "#f87171";
    });

    st.textContent = "STARTING CAMERA + TRACKER...";
    try {
      lmkr = await Promise.race([
        trackerPromise,
        new Promise((_, reject) => setTimeout(
          () => reject(new Error("Tracker download timed out")),
          30000
        ))
      ]);
      st.textContent = "HAND CONTROL READY";
    } catch (e) {
      st.textContent = "TRACKER LOAD FAILED - RETRY";
      st.title = e && e.message ? e.message : "Tracker unavailable";
      st.style.color = "#f87171";
      st.style.cursor = "pointer";
      trackerPromise.then(loaded => {
        lmkr = loaded;
        if (stream) {
          st.textContent = "HAND CONTROL READY";
          st.style.color = "#22d3ee";
          st.style.cursor = "";
        }
      }).catch(() => {});
    }
  } catch (e) {
    st.textContent = "GESTURE CONTROL BLOCKED";
    st.title = e && e.message ? e.message : "Unknown camera or model error";
    st.style.color = "#f87171";
  }
}
st.addEventListener("click", async () => {
  if (lmkr || st.textContent !== "TRACKER LOAD FAILED - RETRY") return;
  st.textContent = "RETRYING TRACKER...";
  st.title = "";
  st.style.color = "#22d3ee";
  try {
    lmkr = await makeLandmarker();
    st.textContent = "HAND CONTROL READY";
    st.style.cursor = "";
    if (stream) loop();
  } catch (e) {
    st.textContent = "TRACKER LOAD FAILED - RETRY";
    st.title = e && e.message ? e.message : "Tracker unavailable";
    st.style.color = "#f87171";
    st.style.cursor = "pointer";
  }
});
function loop() {
  const now = performance.now();
  if (lmkr && vid.readyState >= 2 && vid.currentTime !== lastT && now - lastDetect >= 70) {
    lastDetect = now;
    lastT = vid.currentTime;
    let res;
    try {
      res = lmkr.detectForVideo(vid, performance.now());
    } catch (e) {
      st.textContent = "HAND TRACKER PAUSED";
      st.style.color = "#f87171";
      requestAnimationFrame(loop);
      return;
    }
    ctx.clearRect(0, 0, cvs.width, cvs.height);
    const hands = res.landmarks || [];
    ctx.fillStyle = "#22d3ee";
    hands.forEach(h => h.forEach(p => { ctx.beginPath(); ctx.arc(p.x * cvs.width, p.y * cvs.height, 2.5, 0, 6.28); ctx.fill(); }));
    if (hands.length === 2) {
      prev = null;
      const c1 = palmCenter(hands[0]), c2 = palmCenter(hands[1]);
      const rawD = dist(c1, c2);
      smoothZoomDistance = smoothZoomDistance === null ? rawD : smoothZoomDistance * 0.72 + rawD * 0.28;
      if (prevD !== null && Math.abs(smoothZoomDistance - prevD) > 0.012 && window.ctrl3d) {
        window.ctrl3d.zoom(1 + (smoothZoomDistance - prevD) * 2.8);
        st.textContent = smoothZoomDistance > prevD ? "ZOOM IN" : "ZOOM OUT";
      }
      prevD = smoothZoomDistance;
    } else if (hands.length === 1) {
      prevD = null;
      smoothZoomDistance = null;
      const h = hands.at(0), t = h.at(4), idx = h.at(8);
      const indexUp = fingerExtended(h, 5, 6, 8);
      const middleUp = fingerExtended(h, 9, 10, 12);
      const ringUp = fingerExtended(h, 13, 14, 16);
      const pinkyUp = fingerExtended(h, 17, 18, 20);
      const peace = indexUp && middleUp && !ringUp && !pinkyUp;
      if (peace && Date.now() - lastR > 1500) {
        lastR = Date.now(); prev = null;
        if (window.ctrl3d) window.ctrl3d.reset();
        st.textContent = "VIEW RESET";
      } else if (dist(t, idx) < 0.09) {
        const cx = (t.x + idx.x) * 0.5, cy = (t.y + idx.y) * 0.5;
        const rawPalm = {x: cx, y: cy};
        smoothPinch = smoothPinch === null
          ? rawPalm
          : {x: smoothPinch.x * 0.58 + rawPalm.x * 0.42, y: smoothPinch.y * 0.58 + rawPalm.y * 0.42};
        if (prev && window.ctrl3d) {
          const dx = smoothPinch.x - prev.x;
          const dy = smoothPinch.y - prev.y;
          window.ctrl3d.rotate(-dx * 420, dy * 420);
          st.textContent = "ROTATING 3D";
        }
        prev = { x: smoothPinch.x, y: smoothPinch.y };
      } else {
        prev = null;
        smoothPinch = null;
        st.textContent = "PINCH TO ROTATE";
      }
    } else { prev = null; prevD = null; st.textContent = "SHOW HAND"; }
  }
  requestAnimationFrame(loop);
}
window.addEventListener("pagehide", () => {
  if (stream) stream.getTracks().forEach(track => track.stop());
  if (lmkr && typeof lmkr.close === "function") lmkr.close();
});
boot();
</script>"""


def _chem_html(mol: str, port: int | None = None) -> str:
    raw_sdf = _get_molecule_sdf(mol)
    sdf = raw_sdf.replace("\\", "\\\\").replace("`", "\\`")
    title = html.escape(mol)
    return f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <title>{title} - 3D Study Lab</title>
  <style>
    html, body {{ margin:0; width:100%; height:100%; overflow:hidden; background:#020617; }}
    #viewer {{ width:100%; height:100%; }}
  </style>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js"></script>
</head>
<body>
  <div id="viewer"></div>
  <script>
    window.ctrl3d = null;
    const viewer = $3Dmol.createViewer(document.getElementById("viewer"), {{
      backgroundColor: "#020617"
    }});
    viewer.addModel(`{sdf}`, "sdf");
    viewer.setStyle({{}}, {{stick: {{radius:0.12}}, sphere: {{scale:0.25}}}});
    viewer.zoomTo();
    viewer.render();
    window.ctrl3d = {{
      rotate(dx, dy) {{
        viewer.rotate(dx, {{x:0, y:1, z:0}});
        viewer.rotate(dy, {{x:1, y:0, z:0}});
        viewer.render();
      }},
      zoom(factor) {{
        viewer.zoom(factor);
        viewer.render();
      }},
      reset() {{
        viewer.zoomTo();
        viewer.render();
      }}
    }};
    let pointer = null;
    viewer.getCanvas().addEventListener("pointerdown", e => {{
      pointer = {{x:e.clientX, y:e.clientY}};
      viewer.getCanvas().setPointerCapture(e.pointerId);
    }});
    viewer.getCanvas().addEventListener("pointermove", e => {{
      if (!pointer) return;
      const dx = e.clientX - pointer.x;
      const dy = e.clientY - pointer.y;
      window.ctrl3d.rotate(dx * 0.35, dy * 0.35);
      pointer = {{x:e.clientX, y:e.clientY}};
    }});
    viewer.getCanvas().addEventListener("pointerup", () => pointer = null);
    viewer.getCanvas().addEventListener("wheel", e => {{
      e.preventDefault();
      window.ctrl3d.zoom(e.deltaY < 0 ? 1.08 : 0.92);
    }}, {{passive:false}});
  </script>
  {_hud_script(port)}
</body>
</html>"""


def _physics_html(topic: str, port: int | None = None) -> str:
    key = re.sub(r"[^a-z0-9]+", "_", topic.lower()).strip("_")
    title = html.escape(topic)
    descriptions = {
        "wave_interference": "Two coherent waves overlap. Observe constructive and destructive interference.",
        "magnetic_helix": "A magnetic dipole field rendered as rotating helical field lines.",
        "electric_field": "Electric field vectors around a positive and negative charge.",
    }
    description = descriptions.get(key, "Interactive physics visualization.")
    return f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <title>{title} - 3D Study Lab</title>
  <style>
    html, body {{ margin:0; width:100%; height:100%; overflow:hidden; background:#020617; }}
    canvas {{ width:100%; height:100%; display:block; }}
    #label {{
      position:fixed; left:16px; top:16px; max-width:min(420px,70vw);
      color:#e2e8f0; background:rgba(15,23,42,.88); border:1px solid #06b6d4;
      border-radius:10px; padding:10px 14px; font-family:sans-serif; z-index:20;
    }}
  </style>
</head>
<body>
  <div id="label"><strong>{title}</strong><br>{html.escape(description)}</div>
  <canvas id="scene"></canvas>
  {_hud_script(port)}
  <script>
    const canvas = document.getElementById("scene");
    const ctx = canvas.getContext("2d");
    const mode = {json_text(key)};
    let angle = 0, zoom = 1;
    let targetX = 0, targetY = 0;

    function resize() {{
      const dpr = window.devicePixelRatio || 1;
      canvas.width = innerWidth * dpr;
      canvas.height = innerHeight * dpr;
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    }}
    addEventListener("resize", resize);
    resize();

    window.ctrl3d = {{
      rotate(dx, dy) {{
        targetX += dx * 0.01;
        targetY = Math.max(-1.2, Math.min(1.2, targetY + dy * 0.01));
      }},
      zoom(factor) {{
        zoom = Math.max(0.55, Math.min(2.2, zoom * factor));
      }},
      reset() {{
        angle = 0; zoom = 1; targetX = 0; targetY = 0;
      }}
    }};
    let pointer = null;
    canvas.addEventListener("pointerdown", e => {{
      pointer = {{x:e.clientX, y:e.clientY}};
      canvas.setPointerCapture(e.pointerId);
    }});
    canvas.addEventListener("pointermove", e => {{
      if (!pointer) return;
      window.ctrl3d.rotate(
        (e.clientX - pointer.x) * 0.8,
        (e.clientY - pointer.y) * 0.8
      );
      pointer = {{x:e.clientX, y:e.clientY}};
    }});
    canvas.addEventListener("pointerup", () => pointer = null);
    canvas.addEventListener("wheel", e => {{
      e.preventDefault();
      window.ctrl3d.zoom(e.deltaY < 0 ? 1.08 : 0.92);
    }}, {{passive:false}});

    function waveInterference(w, h) {{
      // Apply the same gesture-driven transform used by the other physics scenes.
      ctx.save();
      ctx.translate(w / 2, h / 2);
      ctx.scale(zoom, zoom);
      ctx.rotate(targetX);
      ctx.transform(1, 0, Math.sin(targetY), Math.cos(targetY), 0, 0);
      const halfW = w / 2, halfH = h / 2;
      for (let y = -halfH + 80; y < halfH - 80; y += 12) {{
        for (let x = -halfW + 80; x < halfW - 80; x += 12) {{
          const r1 = Math.hypot(x + halfW * 0.14, y);
          const r2 = Math.hypot(x - halfW * 0.14, y);
          const v = Math.sin(r1 * 0.045 - angle) + Math.sin(r2 * 0.045 - angle);
          const a = Math.min(1, Math.abs(v) * 0.38);
          ctx.fillStyle = v >= 0 ? `rgba(34,211,238,${{a}})` : `rgba(244,114,182,${{a}})`;
          ctx.fillRect(x, y, 10, 10);
        }}
      }}
      ctx.strokeStyle = "#22d3ee";
      ctx.beginPath();
      for (let x = -halfW + 80; x < halfW - 80; x += 4) {{
        const y = Math.sin(x * 0.025 - angle) * 55;
        x === -halfW + 80 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
      }}
      ctx.stroke();
      ctx.restore();
    }}

    function magneticHelix(w, h) {{
      ctx.save();
      ctx.translate(w / 2, h / 2);
      ctx.scale(zoom, zoom);
      ctx.rotate(targetX);
      ctx.transform(1, 0, Math.sin(targetY), Math.cos(targetY), 0, 0);
      for (let ring = -280; ring <= 280; ring += 16) {{
        ctx.strokeStyle = `hsla(${{190 + ring / 5}},90%,65%,0.8)`;
        ctx.beginPath();
        for (let s = 0; s <= 360; s += 6) {{
          const t = s * Math.PI / 180;
          const x = ring + 70 * Math.sin(t);
          const y = 70 * Math.cos(t) + ring * 0.18;
          s === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
        }}
        ctx.stroke();
      }}
      ctx.restore();
    }}

    function electricField(w, h) {{
      ctx.save();
      ctx.translate(w / 2, h / 2);
      ctx.scale(zoom, zoom);
      ctx.rotate(targetX);
      ctx.transform(1, 0, Math.sin(targetY), Math.cos(targetY), 0, 0);
      const p1 = {{x:-180,y:0,charge:1}}, p2 = {{x:180,y:0,charge:-1}};
      for (let y = -260; y <= 260; y += 28) {{
        for (let x = -420; x <= 420; x += 28) {{
          let ex = 0, ey = 0;
          for (const p of [p1,p2]) {{
            const dx=x-p.x, dy=y-p.y, r2=dx*dx+dy*dy+1200;
            ex += p.charge*dx/r2;
            ey += p.charge*dy/r2;
          }}
          const len = Math.hypot(ex,ey)||1;
          ex /= len; ey /= len;
          ctx.strokeStyle="rgba(34,211,238,.65)";
          ctx.beginPath(); ctx.moveTo(x-ex*8,y-ey*8); ctx.lineTo(x+ex*8,y+ey*8); ctx.stroke();
        }}
      }}
      ctx.fillStyle="#ef4444"; ctx.beginPath(); ctx.arc(p1.x,p1.y,18,0,Math.PI*2); ctx.fill();
      ctx.fillStyle="#3b82f6"; ctx.beginPath(); ctx.arc(p2.x,p2.y,18,0,Math.PI*2); ctx.fill();
      ctx.restore();
    }}

    function draw() {{
      ctx.fillStyle="#020617"; ctx.fillRect(0,0,innerWidth,innerHeight);
      angle += 0.025;
      if (mode === "magnetic_helix") magneticHelix(innerWidth, innerHeight);
      else if (mode === "electric_field") electricField(innerWidth, innerHeight);
      else waveInterference(innerWidth, innerHeight);
      requestAnimationFrame(draw);
    }}
    draw();
  </script>
</body>
</html>"""


def json_text(value: str) -> str:
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def run(parameters: dict, player=None, session_memory=None) -> str:
    subject = str(parameters.get("subject") or "").strip().lower()
    topic = str(parameters.get("topic") or "").strip()
    if not subject or not topic:
        return "Sir, I need both a subject and a topic for the 3D study model."

    if subject not in ("chemistry", "physics"):
        return "Sir, the 3D study model supports chemistry or physics."

    try:
        key = re.sub(r"[^a-z0-9]+", "_", topic.lower()).strip("_") or "model"
        port = _start_server()
        # Download the hand model in the background; do not delay opening the 3D viewer.
        threading.Thread(target=_prepare_tracker_model, kwargs={'timeout': 20}, daemon=True).start()
        if subject == "chemistry":
            html_document = _chem_html(topic, port)
            filename = f"chemistry_{key}.html"
        else:
            html_document = _physics_html(topic, port)
            filename = f"physics_{key}.html"

        _DIR.mkdir(parents=True, exist_ok=True)
        _CACHE_DIR.mkdir(parents=True, exist_ok=True)
        html_file = _DIR / filename
        html_file.write_text(html_document, encoding="utf-8")

        url = f"http://127.0.0.1:{port}/{urllib.parse.quote(filename)}"
        opened = webbrowser.open(url)
        if player:
            try:
                player.write_log(f"JARVIS: Opened 3D study model for {topic}.")
            except Exception:
                pass

        if not opened:
            return f"Sir, I prepared the 3D study model at {url}, but my browser could not open it."
        return f"Opened the 3D {subject} study model for {topic}."
    except Exception as exc:
        return f"Sir, I could not open the 3D study model: {exc}"
