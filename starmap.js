import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { CSS2DRenderer, CSS2DObject } from 'three/addons/renderers/CSS2DRenderer.js';

// ---------------------------------------------------------------------------
// Config. Coordinates are in parsecs (HYG catalog convention, heliocentric).
// 1 scene unit = 1 parsec. Nearest 1000 stars all fit within ~15.5 pc of Sol.
// ---------------------------------------------------------------------------
const CONFIG = {
  catalog: './stars.json',
  constellations: './constellations.json',
  gridPlaneSize: 2000,        // procedural-grid plane — much larger than camera range
  gridFineSpacing: 2,         // pc — fine grid cells, only visible up close
  gridCoarseSpacing: 20,      // pc — coarse grid cells, visible at all ranges
  gridFineFadeStart: 30,      // pc from camera — fine grid begins fading
  gridFineFadeEnd: 100,       // pc from camera — fine grid fully gone
  gridFadeStart: 60,          // pc from camera — coarse grid begins fading
  gridFadeEnd: 600,           // pc from camera — coarse grid fully gone (so plane edge is unreachable)
  axisHalfLength: 500,
  fogNear: 30,
  fogFar: 400,
  labelFadeNear: 30,
  labelFadeFar: 80,
  maxConnectionDist: 2.0,
  minConnectionDist: 0.05,
  stemBaseOpacity: 0.10,      // per-stem alpha at the reference star count
  stemRefCount: 1500,         // bumped from 1000 so stems fade 1.5x slower as N grows
  bgStarCount: 3500,
  bgRadius: 5000,
  accentColor: 0x00d4ff,
  solColor: 0xffd070,
  starSizeRange: [0.67, 3.11],
  solSize: 3.53,
  magBrightRef: -2,
  magFaintRef: 14,
};

// Standard-ish color for each Morgan-Keenan spectral class.
const SPECT_COLORS = {
  O: 0x9bb0ff, B: 0xaabfff, A: 0xcad7ff, F: 0xfff4ea,
  G: 0xffeaa1, K: 0xffc878, M: 0xff9966, L: 0xff7a55,
  T: 0xc73d3d, Y: 0xa02828, W: 0xb8e8ff, C: 0xff5050,
  '?': 0xb8e8ff,
};

// ---------------------------------------------------------------------------
// Naming heuristic: anything that doesn't look like a Bayer designation,
// Flamsteed number, or catalog id is treated as a "proper" name worth labeling.
// HYG packs proper names (e.g. "Sirius"), Bayer/Flamsteed (e.g. "Alp Cen",
// "104 Tau"), and bare catalog ids ("HIP 12345") all into the same field.
// ---------------------------------------------------------------------------
function isProperName(name) {
  if (!name) return false;
  if (name === 'Sol') return true;
  if (/^(HIP|HD|Gl|GJ|ID|TYC|2MASS)\s/i.test(name)) return false;
  if (/\d/.test(name)) return false; // Flamsteed has numbers; proper names typically don't
  const parts = name.split(/\s+/);
  // 3-letter capitalized constellation abbreviations (e.g. "Alp Cen", "p Eri")
  if (parts.length >= 2) {
    const last = parts[parts.length - 1];
    if (last.length === 3 && /^[A-Z][a-z]{2}$/.test(last)) return false;
  }
  return true;
}

// Apparent magnitude → screen point size. Lower mag = brighter = bigger.
function magToSize(mag) {
  const [maxSize, minSize] = [CONFIG.starSizeRange[1], CONFIG.starSizeRange[0]];
  if (mag == null) return (maxSize + minSize) / 2;
  // Clamp inside the reference window then linearly interpolate.
  const clamped = Math.max(CONFIG.magBrightRef, Math.min(CONFIG.magFaintRef, mag));
  const t = (clamped - CONFIG.magBrightRef) / (CONFIG.magFaintRef - CONFIG.magBrightRef);
  return maxSize + (minSize - maxSize) * t;
}

function colorForSpect(cls) {
  return new THREE.Color(SPECT_COLORS[cls] ?? SPECT_COLORS['?']);
}

// ---------------------------------------------------------------------------
// Build a "soft dot" sprite texture for stars — additive blend gives glow.
// ---------------------------------------------------------------------------
function buildStarTexture() {
  const size = 128;
  const canvas = document.createElement('canvas');
  canvas.width = canvas.height = size;
  const ctx = canvas.getContext('2d');
  const grad = ctx.createRadialGradient(size / 2, size / 2, 0, size / 2, size / 2, size / 2);
  grad.addColorStop(0.00, 'rgba(255,255,255,1.0)');
  grad.addColorStop(0.18, 'rgba(255,255,255,0.85)');
  grad.addColorStop(0.45, 'rgba(255,255,255,0.30)');
  grad.addColorStop(1.00, 'rgba(255,255,255,0.0)');
  ctx.fillStyle = grad;
  ctx.fillRect(0, 0, size, size);
  const tex = new THREE.CanvasTexture(canvas);
  tex.colorSpace = THREE.SRGBColorSpace;
  return tex;
}

// ---------------------------------------------------------------------------
// Each star → its single nearest neighbor (within max distance). Stars are
// already sorted by distance from Sol; we return line verts + a parallel
// array of per-pair "max endpoint distance from Sol", sorted ascending.
// That lets the slider use setDrawRange to show only connections whose far
// endpoint is within the current radius — no per-frame recomputation.
//
// Uses a uniform spatial grid (cell size = maxConnectionDist) so each star
// only checks its 27 neighboring cells. O(n) overall — handles 10k+ stars in
// a few milliseconds instead of seconds.
// ---------------------------------------------------------------------------
function buildConnections(positions /* Float32Array length n*3 */, distances /* Float32Array length n */) {
  const n = positions.length / 3;
  const cellSize = CONFIG.maxConnectionDist;
  const maxD2 = cellSize ** 2;
  const minD2 = CONFIG.minConnectionDist ** 2;

  // Bucket stars by 3D cell.
  const buckets = new Map();
  for (let i = 0; i < n; i++) {
    const cx = Math.floor(positions[i * 3]     / cellSize);
    const cy = Math.floor(positions[i * 3 + 1] / cellSize);
    const cz = Math.floor(positions[i * 3 + 2] / cellSize);
    const key = `${cx},${cy},${cz}`;
    let bucket = buckets.get(key);
    if (!bucket) { bucket = []; buckets.set(key, bucket); }
    bucket.push(i);
  }

  const seen = new Set();
  const pairs = []; // [a, b, maxDistFromSol]
  for (let i = 0; i < n; i++) {
    const ix = positions[i * 3], iy = positions[i * 3 + 1], iz = positions[i * 3 + 2];
    const cx = Math.floor(ix / cellSize);
    const cy = Math.floor(iy / cellSize);
    const cz = Math.floor(iz / cellSize);
    let bestJ = -1, bestD2 = Infinity;
    for (let dx = -1; dx <= 1; dx++) {
      for (let dy = -1; dy <= 1; dy++) {
        for (let dz = -1; dz <= 1; dz++) {
          const bucket = buckets.get(`${cx + dx},${cy + dy},${cz + dz}`);
          if (!bucket) continue;
          for (let b = 0; b < bucket.length; b++) {
            const j = bucket[b];
            if (j === i) continue;
            const ddx = positions[j * 3]     - ix;
            const ddy = positions[j * 3 + 1] - iy;
            const ddz = positions[j * 3 + 2] - iz;
            const d2 = ddx * ddx + ddy * ddy + ddz * ddz;
            if (d2 < minD2) continue;
            if (d2 < bestD2) { bestD2 = d2; bestJ = j; }
          }
        }
      }
    }
    if (bestJ >= 0 && bestD2 <= maxD2) {
      const a = Math.min(i, bestJ), b = Math.max(i, bestJ);
      const key = `${a},${b}`;
      if (seen.has(key)) continue;
      seen.add(key);
      pairs.push([a, b, Math.max(distances[a], distances[b])]);
    }
  }
  pairs.sort((p, q) => p[2] - q[2]);

  const verts = new Float32Array(pairs.length * 6);
  const pairMaxDist = new Float32Array(pairs.length);
  for (let k = 0; k < pairs.length; k++) {
    const [a, b, d] = pairs[k];
    verts[k * 6 + 0] = positions[a * 3];     verts[k * 6 + 1] = positions[a * 3 + 1]; verts[k * 6 + 2] = positions[a * 3 + 2];
    verts[k * 6 + 3] = positions[b * 3];     verts[k * 6 + 4] = positions[b * 3 + 1]; verts[k * 6 + 5] = positions[b * 3 + 2];
    pairMaxDist[k] = d;
  }
  return { verts, pairMaxDist };
}

// ---------------------------------------------------------------------------
// Scene setup.
// ---------------------------------------------------------------------------
const container = document.getElementById('scene');

const scene = new THREE.Scene();
scene.background = new THREE.Color(0x000000);
// Fog fades distant grid lines + axis so they feel like they extend to infinity
// rather than terminating at a visible edge.
// Exponential fog. Density tuned so the cluster stays mostly visible even
// when fully zoomed out (camera up to 500 pc from origin).
scene.fog = new THREE.FogExp2(0x000000, 0.0012);

const camera = new THREE.PerspectiveCamera(55, window.innerWidth / window.innerHeight, 0.05, 20000);
camera.position.set(14, 9, 18);

const renderer = new THREE.WebGLRenderer({ antialias: true, powerPreference: 'high-performance' });
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
renderer.setSize(window.innerWidth, window.innerHeight);
container.appendChild(renderer.domElement);

const labelRenderer = new CSS2DRenderer();
labelRenderer.setSize(window.innerWidth, window.innerHeight);
labelRenderer.domElement.style.position = 'absolute';
labelRenderer.domElement.style.top = '0';
labelRenderer.domElement.style.left = '0';
labelRenderer.domElement.style.pointerEvents = 'none';
container.appendChild(labelRenderer.domElement);

const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
controls.dampingFactor = 0.06;
controls.rotateSpeed = 0.6;
controls.zoomSpeed = 0.8;
controls.enablePan = false;
controls.minDistance = 2;
controls.maxDistance = 500;
controls.minPolarAngle = Math.PI * 0.05;
controls.maxPolarAngle = Math.PI * 0.95;
controls.autoRotate = true;
controls.autoRotateSpeed = 0.117;
controls.target.set(0, 0, 0);

// Current view mode — declared here (above tick()) because tick() reads it
// every frame; the view-toggle wiring further down only mutates it.
let mode = 'ext';

// Pause auto-orbit on interaction; resume after idle.
let autoResumeTimer = null;
controls.addEventListener('start', () => {
  controls.autoRotate = false;
  if (autoResumeTimer) clearTimeout(autoResumeTimer);
});
controls.addEventListener('end', () => {
  autoResumeTimer = setTimeout(() => { controls.autoRotate = true; }, 10500);
});

// ---------------------------------------------------------------------------
// Static scene elements: ecliptic grid, Y axis pin, Sol marker, background.
// ---------------------------------------------------------------------------

// Procedural infinite grid — a single large plane whose fragment shader draws
// grid lines using screen-space derivatives, so lines stay crisp at any zoom
// and fade smoothly with distance. The plane is huge but never reaches its
// edge visually because the shader fade hits zero long before that.
const gridMesh = (() => {
  const geom = new THREE.PlaneGeometry(CONFIG.gridPlaneSize, CONFIG.gridPlaneSize);
  geom.rotateX(-Math.PI / 2); // lay flat in XZ plane

  const mat = new THREE.ShaderMaterial({
    uniforms: {
      uCameraPos: { value: new THREE.Vector3() },
      uColor: { value: new THREE.Color(CONFIG.accentColor) },
      uFineSpacing: { value: CONFIG.gridFineSpacing },
      uCoarseSpacing: { value: CONFIG.gridCoarseSpacing },
      uFineFadeStart: { value: CONFIG.gridFineFadeStart },
      uFineFadeEnd: { value: CONFIG.gridFineFadeEnd },
      uFadeStart: { value: CONFIG.gridFadeStart },
      uFadeEnd: { value: CONFIG.gridFadeEnd },
    },
    vertexShader: `
      varying vec3 vWorldPos;
      void main() {
        vec4 wp = modelMatrix * vec4(position, 1.0);
        vWorldPos = wp.xyz;
        gl_Position = projectionMatrix * viewMatrix * wp;
      }
    `,
    fragmentShader: `
      precision highp float;
      varying vec3 vWorldPos;
      uniform vec3 uCameraPos;
      uniform vec3 uColor;
      uniform float uFineSpacing;
      uniform float uCoarseSpacing;
      uniform float uFineFadeStart;
      uniform float uFineFadeEnd;
      uniform float uFadeStart;
      uniform float uFadeEnd;

      // Screen-space-derivative grid: 1 on the line, 0 between, anti-aliased.
      float gridFactor(vec2 p, float spacing) {
        vec2 coord = p / spacing;
        vec2 grid = abs(fract(coord - 0.5) - 0.5) / fwidth(coord);
        return 1.0 - min(min(grid.x, grid.y), 1.0);
      }

      void main() {
        vec2 p = vWorldPos.xz;
        float fine = gridFactor(p, uFineSpacing);
        float coarse = gridFactor(p, uCoarseSpacing);

        float d = distance(vWorldPos, uCameraPos);
        float fineMask = 1.0 - smoothstep(uFineFadeStart, uFineFadeEnd, d);
        float globalFade = 1.0 - smoothstep(uFadeStart, uFadeEnd, d);

        float alpha = max(fine * fineMask * 0.12, coarse * 0.20) * globalFade;
        if (alpha < 0.003) discard;
        gl_FragColor = vec4(uColor, alpha);
      }
    `,
    transparent: true,
    depthWrite: false,
    fog: false, // we handle distance fade ourselves
    side: THREE.DoubleSide,
  });

  const mesh = new THREE.Mesh(geom, mat);
  scene.add(mesh);
  return mesh;
})();
const axisLine = (() => {
  const h = CONFIG.axisHalfLength;
  const geom = new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(0, -h, 0), new THREE.Vector3(0, h, 0)]);
  const mat = new THREE.LineBasicMaterial({ color: CONFIG.accentColor, transparent: true, opacity: 0.187, depthWrite: false });
  const line = new THREE.Line(geom, mat);
  scene.add(line);
  return line;
})();

// Background "infinity" starfield. Parented logically to the camera each
// frame (we copy camera.position into it) so orbit movement has zero parallax.
const bgStars = (() => {
  const positions = new Float32Array(CONFIG.bgStarCount * 3);
  // Stable PRNG (mulberry32) so background layout is reproducible.
  let s = 0xC0FFEE >>> 0;
  const rand = () => {
    s = (s + 0x6D2B79F5) >>> 0;
    let t = Math.imul(s ^ (s >>> 15), 1 | s);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
  for (let i = 0; i < CONFIG.bgStarCount; i++) {
    const theta = 2 * Math.PI * rand();
    const phi = Math.acos(2 * rand() - 1);
    positions[i * 3 + 0] = CONFIG.bgRadius * Math.sin(phi) * Math.cos(theta);
    positions[i * 3 + 1] = CONFIG.bgRadius * Math.sin(phi) * Math.sin(theta);
    positions[i * 3 + 2] = CONFIG.bgRadius * Math.cos(phi);
  }
  const geom = new THREE.BufferGeometry();
  geom.setAttribute('position', new THREE.BufferAttribute(positions, 3));
  const mat = new THREE.PointsMaterial({
    color: 0xffffff, size: 1.4, sizeAttenuation: false,
    transparent: true, opacity: 0.18, depthWrite: false,
    fog: false, // background starfield is "at infinity" — must not get fogged out
  });
  const points = new THREE.Points(geom, mat);
  scene.add(points);
  return points;
})();

// ---------------------------------------------------------------------------
// Constellations. Pre-computed line endpoints from a build-time HYG cross-
// reference (see scripts that produced constellations.json). The endpoints
// are world-space coords just like the star positions; we render them as a
// single LineSegments mesh with a slightly brighter alpha than the regular
// nearest-neighbor connection lines. Visible in both view modes.
// ---------------------------------------------------------------------------
async function loadConstellations() {
  try {
    const res = await fetch(CONFIG.constellations);
    if (!res.ok) return null;
    return await res.json();
  } catch { return null; }
}

let constellationCtx = null;

function buildConstellationLines(data) {
  if (!data || !data.figures) return null;

  // Two parallel things to build:
  //   1. A single LineSegments mesh of all figure segments, sorted by farther-
  //      endpoint distance so the slider can use setDrawRange to hide lines
  //      whose endpoints aren't currently rendered as stars.
  //   2. A CSS2D label per figure positioned at its centroid, so each figure
  //      is identified by name (Boötes &c.) on the map.
  const items = []; // each: { s, max }
  const labels = []; // each: { obj, maxDist } — for slider filtering + toggle
  for (const fig of data.figures) {
    let cx = 0, cy = 0, cz = 0, count = 0, figureMax = 0;
    for (const s of fig.segments) {
      const d1 = Math.hypot(s[0], s[1], s[2]);
      const d2 = Math.hypot(s[3], s[4], s[5]);
      items.push({ s, max: Math.max(d1, d2) });
      cx += s[0] + s[3]; cy += s[1] + s[4]; cz += s[2] + s[5];
      count += 2;
      figureMax = Math.max(figureMax, d1, d2);
    }
    if (!count) continue;
    cx /= count; cy /= count; cz /= count;
    const el = document.createElement('div');
    el.className = 'constellation-label';
    el.textContent = fig.name;
    const obj = new CSS2DObject(el);
    obj.position.set(cx, cy, cz);
    scene.add(obj);
    // The figure's label shows once the slider radius reaches the figure's
    // farthest endpoint — i.e. when the full figure is drawn.
    labels.push({ obj, maxDist: figureMax });
  }
  if (!items.length) return null;
  items.sort((a, b) => a.max - b.max);

  const verts = new Float32Array(items.length * 6);
  const segMaxDist = new Float32Array(items.length);
  for (let i = 0; i < items.length; i++) {
    const s = items[i].s;
    verts[i * 6 + 0] = s[0]; verts[i * 6 + 1] = s[1]; verts[i * 6 + 2] = s[2];
    verts[i * 6 + 3] = s[3]; verts[i * 6 + 4] = s[4]; verts[i * 6 + 5] = s[5];
    segMaxDist[i] = items[i].max;
  }

  const geom = new THREE.BufferGeometry();
  geom.setAttribute('position', new THREE.BufferAttribute(verts, 3));
  const mat = new THREE.LineBasicMaterial({
    color: CONFIG.accentColor,
    transparent: true,
    opacity: 0.55,
    depthWrite: false,
    blending: THREE.AdditiveBlending,
  });
  const mesh = new THREE.LineSegments(geom, mat);
  scene.add(mesh);
  constellationCtx = { mesh, segMaxDist, labels };

  // Apply current slider state immediately so figures don't briefly show in
  // full before the star catalog loads and the first applyRange fires.
  const slider = document.getElementById('range-slider');
  applyConstellationRange(slider ? parseFloat(slider.value) : 0);

  return mesh;
}

function applyConstellationRange(radius) {
  if (!constellationCtx) return;
  const { mesh, segMaxDist, labels } = constellationCtx;
  let cCount = 0;
  while (cCount < segMaxDist.length && segMaxDist[cCount] <= radius) cCount++;
  mesh.geometry.setDrawRange(0, cCount * 2);
  const showLines = mesh.visible;
  for (const { obj, maxDist } of labels) {
    obj.visible = showLines && maxDist <= radius;
  }
}

// ---------------------------------------------------------------------------
// Load catalog and build star geometries.
//
// loadJSONStreaming pulls the response body chunk-by-chunk so we can drive
// a progress bar. If the server emits Content-Encoding: gzip|br, browsers
// hand back decompressed bytes through the reader while Content-Length is
// the COMPRESSED total — we scale the expected total by a typical JSON
// compression ratio so the bar still tracks something useful, and the
// caller clamps the displayed fraction anyway.
// ---------------------------------------------------------------------------
async function loadJSONStreaming(url, onProgress) {
  const res = await fetch(url);
  if (!res.ok) throw new Error(`Failed to load ${url}: ${res.status}`);
  if (!res.body || !res.body.getReader) return res.json();

  const headerTotal = +res.headers.get('content-length') || 0;
  const enc = (res.headers.get('content-encoding') || '').toLowerCase();
  const estimatedTotal = headerTotal * (enc === 'gzip' ? 3 : enc === 'br' ? 4 : 1);

  const reader = res.body.getReader();
  const chunks = [];
  let received = 0;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    chunks.push(value);
    received += value.length;
    if (onProgress) onProgress(received, estimatedTotal);
  }
  // Final progress tick with actual size so the bar settles at 100%.
  if (onProgress) onProgress(received, received);

  const buf = new Uint8Array(received);
  let off = 0;
  for (const c of chunks) { buf.set(c, off); off += c.length; }
  return JSON.parse(new TextDecoder().decode(buf));
}

function loadCatalog(onProgress) {
  return loadJSONStreaming(CONFIG.catalog, onProgress);
}

function buildStars(payload) {
  const rawStars = payload.stars;
  const n = rawStars.length;

  // Sort by distance from Sol so we can use setDrawRange to "scrub" stars
  // in/out as the slider moves. Sol is at the origin so it sorts to index 0.
  const indexed = rawStars.map((row) => {
    const x = row[1], y = row[2], z = row[3];
    return { row, dist: Math.sqrt(x * x + y * y + z * z) };
  });
  indexed.sort((a, b) => a.dist - b.dist);

  const positions = new Float32Array(n * 3);
  const colors = new Float32Array(n * 3);
  const sizes = new Float32Array(n);
  const stemVerts = new Float32Array(n * 6);
  const distances = new Float32Array(n);

  const solCol = new THREE.Color(CONFIG.solColor);

  for (let i = 0; i < n; i++) {
    const { row, dist } = indexed[i];
    const [name, x, y, z, mag, spect] = row;
    distances[i] = dist;
    positions[i * 3] = x; positions[i * 3 + 1] = y; positions[i * 3 + 2] = z;
    const col = name === 'Sol' ? solCol : colorForSpect(spect);
    colors[i * 3] = col.r; colors[i * 3 + 1] = col.g; colors[i * 3 + 2] = col.b;
    sizes[i] = name === 'Sol' ? CONFIG.solSize : magToSize(mag);
    stemVerts[i * 6] = x;     stemVerts[i * 6 + 1] = y; stemVerts[i * 6 + 2] = z;
    stemVerts[i * 6 + 3] = x; stemVerts[i * 6 + 4] = 0; stemVerts[i * 6 + 5] = z;
  }

  const tex = buildStarTexture();

  // Stars — single Points cloud with custom shader (per-vertex size + color).
  const starGeom = new THREE.BufferGeometry();
  starGeom.setAttribute('position', new THREE.BufferAttribute(positions, 3));
  starGeom.setAttribute('color', new THREE.BufferAttribute(colors, 3));
  starGeom.setAttribute('size', new THREE.BufferAttribute(sizes, 1));
  const starMat = new THREE.ShaderMaterial({
    uniforms: { map: { value: tex } },
    vertexShader: `
      attribute float size;
      attribute vec3 color;
      varying vec3 vColor;
      // Cap the absolute pixel size — without this, stars within a few pc of
      // the camera (Sirius / α Cen / etc. in FROM-SOL view) explode into
      // huge blobs because of the 1/z scaling. Distant stars still shrink
      // normally, so the far sky stays satisfyingly full.
      const float MAX_PX = 14.0;
      void main() {
        vColor = color;
        vec4 mv = modelViewMatrix * vec4(position, 1.0);
        gl_PointSize = min(MAX_PX, size * (300.0 / -mv.z));
        gl_Position = projectionMatrix * mv;
      }
    `,
    fragmentShader: `
      uniform sampler2D map;
      varying vec3 vColor;
      void main() {
        vec4 tex = texture2D(map, gl_PointCoord);
        if (tex.a < 0.02) discard;
        gl_FragColor = vec4(vColor, 1.0) * tex;
      }
    `,
    transparent: true, depthWrite: false, blending: THREE.AdditiveBlending,
  });
  const starsPoints = new THREE.Points(starGeom, starMat);
  scene.add(starsPoints);

  // Stems to the ecliptic plane.
  const stemGeom = new THREE.BufferGeometry();
  stemGeom.setAttribute('position', new THREE.BufferAttribute(stemVerts, 3));
  const stemLines = new THREE.LineSegments(stemGeom, new THREE.LineBasicMaterial({
    color: CONFIG.accentColor, transparent: true, opacity: 0.10, depthWrite: false,
  }));
  scene.add(stemLines);

  // Connection lines — pre-sorted by max-endpoint-distance from Sol.
  const { verts: connVerts, pairMaxDist } = buildConnections(positions, distances);
  const connGeom = new THREE.BufferGeometry();
  connGeom.setAttribute('position', new THREE.BufferAttribute(connVerts, 3));
  const connLines = new THREE.LineSegments(connGeom, new THREE.LineBasicMaterial({
    color: CONFIG.accentColor, transparent: true, opacity: 0.18,
    depthWrite: false, blending: THREE.AdditiveBlending,
  }));
  scene.add(connLines);

  // Labels — only for proper-named stars to avoid clutter. Track distance so
  // we can toggle visibility against the slider radius.
  for (let i = 0; i < n; i++) {
    const { row, dist } = indexed[i];
    const [name, x, y, z] = row;
    if (!isProperName(name)) continue;
    const el = document.createElement('div');
    el.className = name === 'Sol' ? 'star-label hub' : 'star-label';
    el.textContent = name.toUpperCase();
    const obj = new CSS2DObject(el);
    obj.position.set(x, y, z);
    scene.add(obj);
    labelData.push({ obj, dist });
  }

  return { distances, pairMaxDist, starsPoints, stemLines, connLines, total: n };
}

// Labels are DOM (CSS2D), so fog doesn't apply — fade them manually with
// distance from camera so they don't stay crisp when the rest of the scene
// has faded to black. Each entry: { obj, dist, lastOpacity }.
const labelData = [];
function updateLabelFade() {
  const near = CONFIG.labelFadeNear;
  const far = CONFIG.labelFadeFar;
  const span = far - near;
  for (const entry of labelData) {
    if (!entry.obj.visible) continue;
    const d = entry.obj.position.distanceTo(camera.position);
    let t = d <= near ? 1 : d >= far ? 0 : 1 - (d - near) / span;
    // Quantize to ~1% so auto-orbit's tiny per-frame deltas don't trigger
    // a DOM write for every label every frame.
    t = Math.round(t * 100) / 100;
    if (entry.lastOpacity !== t) {
      entry.obj.element.style.opacity = t.toFixed(2);
      entry.lastOpacity = t;
    }
  }
}

// Apply the slider radius: rendered star count, stem count, connection count,
// and label visibility are all derived in O(log n) using upper-bound scans
// (linear here for simplicity — only runs on slider input, ~ms for 4k stars).
function applyRange(ctx, radius) {
  let starCount = 0;
  while (starCount < ctx.distances.length && ctx.distances[starCount] <= radius) starCount++;
  let connCount = 0;
  while (connCount < ctx.pairMaxDist.length && ctx.pairMaxDist[connCount] <= radius) connCount++;

  ctx.starsPoints.geometry.setDrawRange(0, starCount);
  ctx.stemLines.geometry.setDrawRange(0, starCount * 2);
  ctx.connLines.geometry.setDrawRange(0, connCount * 2);

  // Stems use normal (non-additive) alpha blending, so dense overlap saturates
  // quickly. Scale per-stem alpha ~ 1/count to hold total visual density near
  // the value tuned at the default count.
  const ratio = CONFIG.stemRefCount / Math.max(starCount, 1);
  ctx.stemLines.material.opacity = CONFIG.stemBaseOpacity * Math.min(1, ratio);

  let labeled = 0;
  for (const entry of labelData) {
    const visible = entry.dist <= radius;
    entry.obj.visible = visible;
    if (visible) labeled++;
  }

  applyConstellationRange(radius);

  const hud = document.getElementById('hud-count');
  if (hud) hud.textContent = `${starCount} SYSTEMS · ${labeled} NAMED`;

  // Invalidate the stationary-camera cache so newly-visible labels get their
  // opacity recomputed on the next tick (even if the camera hasn't moved).
  _prevCamPos.set(Infinity, Infinity, Infinity);
}

// ---------------------------------------------------------------------------
// Resize + render loop.
// ---------------------------------------------------------------------------
function onResize() {
  const w = window.innerWidth;
  const h = window.innerHeight;
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
  renderer.setSize(w, h);
  labelRenderer.setSize(w, h);
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
}
window.addEventListener('resize', onResize);

const _prevCamPos = new THREE.Vector3(Infinity, Infinity, Infinity);
function tick() {
  // OrbitControls.update() drives the camera off its internal spherical state;
  // running it in Sol mode would tug the camera away from the origin where
  // SkyControls just placed it.
  if (mode === 'ext') controls.update();

  // Skip camera-dependent work when the camera didn't actually move this frame
  // (sol view holding still; auto-orbit paused after interaction).
  if (!_prevCamPos.equals(camera.position)) {
    bgStars.position.copy(camera.position);
    gridMesh.material.uniforms.uCameraPos.value.copy(camera.position);
    updateLabelFade();
    _prevCamPos.copy(camera.position);
  }

  renderer.render(scene, camera);
  labelRenderer.render(scene, camera);
  requestAnimationFrame(tick);
}
tick();

// ---------------------------------------------------------------------------
// View modes: ORBIT (default — OrbitControls around Sol) and FROM-SOL (camera
// fixed at the origin, pointer drag rotates yaw/pitch in place). OrbitControls
// can't operate at zero distance from its target so we use a custom controller
// for the Sol view rather than fighting OrbitControls' invariants.
// ---------------------------------------------------------------------------
class SkyControls {
  constructor(camera, dom) {
    this.camera = camera;
    this.dom = dom;
    this.enabled = false;
    this.yaw = 0;
    this.pitch = 0;
    this._down = false;
    this._lx = 0;
    this._ly = 0;
    this._sens = 0.0035;
    this._onDown = (e) => this._handleDown(e);
    this._onMove = (e) => this._handleMove(e);
    this._onUp   = (e) => this._handleUp(e);
  }
  enable() {
    if (this.enabled) return;
    this.enabled = true;
    this.dom.addEventListener('pointerdown', this._onDown);
    this._apply();
  }
  disable() {
    if (!this.enabled) return;
    this.enabled = false;
    this.dom.removeEventListener('pointerdown', this._onDown);
    this.dom.removeEventListener('pointermove', this._onMove);
    this.dom.removeEventListener('pointerup', this._onUp);
    this._down = false;
  }
  _handleDown(e) {
    this._down = true;
    this._lx = e.clientX; this._ly = e.clientY;
    try { this.dom.setPointerCapture(e.pointerId); } catch {}
    this.dom.addEventListener('pointermove', this._onMove);
    this.dom.addEventListener('pointerup', this._onUp);
  }
  _handleMove(e) {
    if (!this._down) return;
    const dx = e.clientX - this._lx;
    const dy = e.clientY - this._ly;
    this._lx = e.clientX; this._ly = e.clientY;
    this.yaw   -= dx * this._sens;
    this.pitch -= dy * this._sens;
    const lim = Math.PI / 2 - 0.01;
    if (this.pitch >  lim) this.pitch =  lim;
    if (this.pitch < -lim) this.pitch = -lim;
    this._apply();
  }
  _handleUp() {
    this._down = false;
    this.dom.removeEventListener('pointermove', this._onMove);
    this.dom.removeEventListener('pointerup', this._onUp);
  }
  _apply() {
    this.camera.position.set(0, 0, 0);
    this.camera.rotation.order = 'YXZ';
    this.camera.rotation.set(this.pitch, this.yaw, 0);
  }
  reset(yaw = 0, pitch = 0) { this.yaw = yaw; this.pitch = pitch; this._apply(); }
}

const skyControls = new SkyControls(camera, renderer.domElement);

const savedExt = { pos: camera.position.clone(), target: controls.target.clone() };

function applyViewMode(newMode) {
  if (newMode === mode) return;

  if (newMode === 'sol') {
    // Save the orbit view so we can restore exactly when toggling back.
    savedExt.pos.copy(camera.position);
    savedExt.target.copy(controls.target);
    controls.enabled = false;
    controls.autoRotate = false;
    if (autoResumeTimer) { clearTimeout(autoResumeTimer); autoResumeTimer = null; }

    gridMesh.visible = false;
    axisLine.visible = false;
    if (sceneCtx) {
      sceneCtx.stemLines.visible = false;
      sceneCtx.connLines.visible = false;
    }

    skyControls.reset(0, 0);
    skyControls.enable();
  } else {
    skyControls.disable();
    camera.position.copy(savedExt.pos);
    controls.target.copy(savedExt.target);
    camera.lookAt(controls.target);
    controls.enabled = true;
    controls.autoRotate = true;

    gridMesh.visible = true;
    axisLine.visible = true;
    if (sceneCtx) {
      sceneCtx.stemLines.visible = true;
      sceneCtx.connLines.visible = true;
    }
  }
  mode = newMode;

  for (const btn of document.querySelectorAll('.hud-toggle-btn')) {
    const active = btn.dataset.view === newMode;
    btn.classList.toggle('active', active);
    btn.setAttribute('aria-selected', String(active));
  }

  const hint = document.getElementById('hud-hint');
  if (hint) hint.textContent = newMode === 'sol'
    ? 'DRAG TO LOOK AROUND'
    : 'DRAG TO ORBIT · PINCH/SCROLL TO ZOOM';
}

for (const btn of document.querySelectorAll('.hud-toggle-btn[data-view]')) {
  btn.addEventListener('click', () => applyViewMode(btn.dataset.view));
}

{
  const cBtn = document.getElementById('toggle-constellations');
  if (cBtn) {
    cBtn.addEventListener('click', () => {
      const active = !cBtn.classList.contains('active');
      cBtn.classList.toggle('active', active);
      cBtn.setAttribute('aria-pressed', String(active));
      if (constellationCtx) {
        constellationCtx.mesh.visible = active;
        const slider = document.getElementById('range-slider');
        applyConstellationRange(slider ? parseFloat(slider.value) : 0);
      }
    });
  }
}

let sceneCtx = null;

// Kick off both loads in parallel so the constellation JSON is fetched while
// the bigger star catalog is downloading — but only *add* the constellations
// to the scene once stars exist, so the iconic figures don't appear hanging
// in empty space during the load window.
const constellationsPromise = loadConstellations();

const loadingEl     = document.getElementById('hud-loading');
const loadingFill   = document.getElementById('hud-loading-fill');
const loadingMeta   = document.getElementById('hud-loading-meta');
function onCatalogProgress(received, total) {
  const frac = total > 0 ? Math.min(1, received / total) : 0;
  if (loadingFill) loadingFill.style.width = (frac * 100).toFixed(1) + '%';
  if (loadingMeta) {
    const r = (received / 1048576).toFixed(1);
    const t = total > 0 ? (total / 1048576).toFixed(1) : '?';
    loadingMeta.textContent = `${r} / ${t} MB`;
  }
}

loadCatalog(onCatalogProgress)
  .then(async (payload) => {
    if (loadingEl) {
      loadingEl.classList.add('is-done');
      setTimeout(() => { loadingEl.style.display = 'none'; }, 400);
    }
    const ctx = buildStars(payload);
    sceneCtx = ctx;
    if (mode === 'sol') {
      ctx.stemLines.visible = false;
      ctx.connLines.visible = false;
    }
    const meta = payload.meta || { r_1000: 15.43, r_max: 30.86 };

    // Configure slider from catalog metadata (so the page works for whatever
    // size catalog we ship — 1000 / 5000 / whatever).
    const slider = document.getElementById('range-slider');
    const valueEl = document.getElementById('range-value');
    slider.min = '0';
    slider.max = String(meta.r_max);
    // ~400 increments across the range — feels smooth on desktop and mobile.
    slider.step = String(Math.max(0.05, meta.r_max / 400));
    slider.value = String(meta.r_1000);

    const update = () => {
      const r = parseFloat(slider.value);
      valueEl.textContent = r.toFixed(1);
      applyRange(ctx, r);
    };
    slider.addEventListener('input', update);
    update();

    // Constellations: build now that stars exist so the figures don't appear
    // hanging in empty space during the catalog download.
    const constellationData = await constellationsPromise;
    if (constellationData) {
      buildConstellationLines(constellationData);
      update(); // re-apply slider to filter the just-built constellation set
    }
  })
  .catch((err) => {
    console.error(err);
    const hudEl = document.getElementById('hud-count');
    if (hudEl) hudEl.textContent = 'CATALOG LOAD FAILED';
  });
