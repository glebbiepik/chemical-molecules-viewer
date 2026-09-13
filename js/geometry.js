// geometry.js — 3D coordinate generation + force-field relaxation
// Strategy: seed coordinates via depth-first placement using ideal bond angles
// (sp3 109.5, sp2 120, sp 180, ring closure fallback), then relax with a
// simple Urey-Bradley-style force field (bonds, angles, nonbonded repulsion).

function vec3(x, y, z) { return { x, y, z }; }
function vsub(a, b) { return vec3(a.x - b.x, a.y - b.y, a.z - b.z); }
function vadd(a, b) { return vec3(a.x + b.x, a.y + b.y, a.z + b.z); }
function vscale(a, s) { return vec3(a.x * s, a.y * s, a.z * s); }
function vlen(a) { return Math.sqrt(a.x * a.x + a.y * a.y + a.z * a.z); }
function vnorm(a) { const l = vlen(a) || 1; return vscale(a, 1 / l); }
function vcross(a, b) {
  return vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
function vdot(a, b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

// Rotate point p around axis (unit) by angle, through origin-relative pivot.
function rotateAroundAxis(p, axis, angle, pivot) {
  const rel = vsub(p, pivot);
  const c = Math.cos(angle), s = Math.sin(angle);
  const a = vnorm(axis);
  // Rodrigues rotation
  const cross = vcross(a, rel);
  const dot = vdot(a, rel);
  const rotated = vadd(
    vadd(vscale(rel, c), vscale(cross, s)),
    vscale(a, dot * (1 - c))
  );
  return vadd(pivot, rotated);
}

function geometryOf(el, bondSum, nHeavy) {
  // returns {angles: [deg], shape: 'linear'|'trig'|'tetra'|'bent'...}
  if (nHeavy <= 1) return 109.47;
  // sp: 2 connections & bond sum suggests triple or double bonds on both sides
  if (nHeavy === 2) {
    if (bondSum >= 3) return 180;     // sp (triple bond, or two doubles: CO2)
    return 120;                        // sp2 (e.g. C=C, C=O)
  }
  if (nHeavy === 3) {
    if (bondSum >= 4) return 120;      // sp2 with one double bond
    return 109.47;                     // sp3 amine / NH3-like
  }
  if (nHeavy === 4) {
    if (bondSum >= 5) return 109.47;   // hypervalent; keep tetrahedral
    return 109.47;                     // sp3 carbon
  }
  return 109.47;
}

// Build initial 3D coordinates for a Molecule.
function build3D(mol) {
  const n = mol.atoms.length;
  const pos = mol.atoms.map(() => vec3(0, 0, 0));
  const placed = new Array(n).fill(false);
  const adj = mol.atoms.map((_, i) => mol.neighbors(i));

  if (n === 0) return pos;

  // BFS from atom 0 (prefer a heavy atom with degree > 1 as the seed)
  let root = 0;
  for (let i = 0; i < n; i++) {
    if (mol.atoms[i].el !== 'H' && adj[i].length > 1) { root = i; break; }
  }
  placed[root] = true;
  pos[root] = vec3(0, 0, 0);

  const queue = [root];
  while (queue.length) {
    const cur = queue.shift();
    const nbs = adj[cur].filter(nb => !placed[nb.j]);
    if (!nbs.length) continue;
    const heavy = nbs.filter(nb => mol.atoms[nb.j].el !== 'H');
    const order = heavy.concat(nbs.filter(nb => mol.atoms[nb.j].el === 'H'));
    const k = order.length;

    // reference direction: away from already-placed neighbors (or arbitrary for root)
    let refDir;
    const placedNbs = adj[cur].filter(nb => placed[nb.j]);
    if (placedNbs.length) {
      refDir = vec3(0, 0, 0);
      for (const nb of placedNbs) refDir = vadd(refDir, vsub(pos[cur], pos[nb.j]));
      refDir = vnorm(refDir);
    } else {
      refDir = vec3(1, 0, 0);
    }

    // perpendicular basis
    let up = Math.abs(refDir.y) < 0.9 ? vec3(0, 1, 0) : vec3(1, 0, 0);
    let t1 = vnorm(vcross(refDir, up));
    let t2 = vnorm(vcross(refDir, t1));

    const deg = adj[cur].length;
    const bs = mol.bondOrderSum(cur);
    const angleDeg = geometryOf(mol.atoms[cur].el, bs, deg);

    for (let m = 0; m < k; m++) {
      const nb = order[m];
      const bl = bondLength(mol.atoms[cur].el, mol.atoms[nb.j].el, nb.order, nb.bond && mol.atoms[nb.j].aromatic && mol.atoms[cur].aromatic);
      let dir;
      if (deg === 1) {
        dir = refDir;
      } else if (deg === 2 && angleDeg === 180) {
        // linear: place opposite existing neighbor
        dir = m === 0 ? refDir : vscale(refDir, -1);
      } else {
        // distribute around cone at ideal angle from -refDir
        // For m children: spread evenly in azimuth
        const theta = (angleDeg * Math.PI) / 180; // angle from refDir
        const phi = (m / Math.max(1, k)) * 2 * Math.PI + (placed[cur] ? 0 : 0.7);
        const sinT = Math.sin(theta), cosT = Math.cos(theta);
        dir = vnorm(vadd(
          vscale(refDir, cosT),
          vadd(vscale(t1, sinT * Math.cos(phi)), vscale(t2, sinT * Math.sin(phi)))
        ));
      }
      pos[nb.j] = vadd(pos[cur], vscale(dir, bl));
      placed[nb.j] = true;
      queue.push(nb.j);
    }
  }

  // Any unplaced (disconnected) atoms: scatter
  for (let i = 0; i < n; i++) {
    if (!placed[i]) {
      pos[i] = vec3(3 + i * 2.2, 0, 0);
      placed[i] = true;
    }
  }
  return pos;
}

// ---------- Force-field relaxation ----------
// Terms: bond stretch (harmonic), angle bend (harmonic via 1-3 distance),
// nonbonded repulsion (soft sphere, ignores 1-2 and 1-3 pairs).
function relax(mol, pos, iterations) {
  const n = mol.atoms.length;
  const iters = iterations || 400;
  const adjSets = mol.atoms.map((_, i) => new Set([i]));
  for (const b of mol.bonds) { adjSets[b.i].add(b.j); adjSets[b.j].add(b.i); }
  // 1-3 pairs
  const pairs13 = new Set();
  for (let i = 0; i < n; i++) {
    for (const nb1 of adjSets[i]) for (const nb2 of adjSets[i]) {
      if (nb1 < nb2) pairs13.add(nb1 * n + nb2);
    }
  }
  const is12 = new Set();
  for (const b of mol.bonds) is12.add(Math.min(b.i, b.j) * n + Math.max(b.i, b.j));

  // Precompute ideal 1-3 distances from ideal angles
  const ideal13 = {};
  for (let i = 0; i < n; i++) {
    const nbs = [...adjSets[i]].filter(j => j !== i);
    const bs = mol.bondOrderSum(i);
    const deg = nbs.length;
    const ang = (geometryOf(mol.atoms[i].el, bs, deg) * Math.PI) / 180;
    for (let a = 0; a < nbs.length; a++) for (let b = a + 1; b < nbs.length; b++) {
      const j = nbs[a], k2 = nbs[b];
      const d1 = bondLength(mol.atoms[i].el, mol.atoms[j].el, mol.bondBetween(i, j) ? mol.bondBetween(i, j).order : 1);
      const d2 = bondLength(mol.atoms[i].el, mol.atoms[k2].el, mol.bondBetween(i, k2) ? mol.bondBetween(i, k2).order : 1);
      const d13 = Math.sqrt(d1 * d1 + d2 * d2 - 2 * d1 * d2 * Math.cos(ang));
      ideal13[Math.min(j, k2) * n + Math.max(j, k2)] = d13;
    }
  }

  const kb = 2.0, ka = 0.6, kr = 0.6;
  const step0 = 0.02;
  for (let it = 0; it < iters; it++) {
    const force = mol.atoms.map(() => vec3(0, 0, 0));

    // bonds
    for (const b of mol.bonds) {
      const a1 = mol.atoms[b.i], a2 = mol.atoms[b.j];
      const d = vsub(pos[b.j], pos[b.i]);
      const len = vlen(d) || 1e-6;
      const ideal = bondLength(a1.el, a2.el, b.order, a1.aromatic && a2.aromatic);
      const f = kb * (len - ideal) * 0.5;
      const dir = vscale(d, 1 / len);
      force[b.i] = vadd(force[b.i], vscale(dir, f));
      force[b.j] = vsub(force[b.j], vscale(dir, f));
    }

    // angles (via 1-3 distances)
    for (const key in ideal13) {
      const ki = +key;
      const j = Math.floor(ki / n), k2 = ki % n;
      const d = vsub(pos[k2], pos[j]);
      const len = vlen(d) || 1e-6;
      const ideal = ideal13[key];
      const f = ka * (len - ideal) * 0.5;
      const dir = vscale(d, 1 / len);
      force[j] = vadd(force[j], vscale(dir, f));
      force[k2] = vsub(force[k2], vscale(dir, f));
    }

    // nonbonded repulsion
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const key = i * n + j;
        if (is12.has(key) || pairs13.has(key)) continue;
        const d = vsub(pos[j], pos[i]);
        const len = vlen(d) || 1e-6;
        const ri = ELEMENTS[mol.atoms[i].el] ? ELEMENTS[mol.atoms[i].el].rcov : 1.2;
        const rj = ELEMENTS[mol.atoms[j].el] ? ELEMENTS[mol.atoms[j].el].rcov : 1.2;
        const rsum = (ri + rj) * 1.3;
        if (len < rsum) {
          const f = kr * (rsum - len) / rsum;
          const dir = vscale(d, 1 / len);
          force[i] = vsub(force[i], vscale(dir, f));
          force[j] = vadd(force[j], vscale(dir, f));
        }
      }
    }

    // integrate with decaying step
    const step = step0 * (1 - it / iters * 0.9) + 0.002;
    for (let i = 0; i < n; i++) {
      pos[i] = vadd(pos[i], vscale(force[i], step));
    }
  }

  // center
  let cx = 0, cy = 0, cz = 0;
  for (const p of pos) { cx += p.x; cy += p.y; cz += p.z; }
  cx /= n; cy /= n; cz /= n;
  for (const p of pos) { p.x -= cx; p.y -= cy; p.z -= cz; }
  return pos;
}

// Full pipeline: SMILES -> Molecule with 3D coords
function moleculeFromSMILES(smiles) {
  const mol = parseSMILES(smiles);
  const pos = build3D(mol);
  relax(mol, pos, 1500);
  mol.pos = pos;
  mol.formula = formulaString(countAtoms(mol));
  mol.molarMass = mol.atoms.reduce((s, a) => s + (ELEMENTS[a.el] ? ELEMENTS[a.el].m : 0), 0);
  return mol;
}

function countAtoms(mol) {
  const c = {};
  for (const a of mol.atoms) c[a.el] = (c[a.el] || 0) + 1;
  return c;
}

// Try to interpret user input as name, formula, or SMILES.
function resolveInput(input) {
  const q = input.trim();
  if (!q) throw new Error('Enter a molecule name, formula, or SMILES');
  const lower = q.toLowerCase().replace(/[\s_-]+/g, ' ').trim();
  // 1. database lookup by name or formula
  for (const key in DB) {
    const [label, smiles] = DB[key];
    const names = label.split(',')[0].trim().toLowerCase();
    const formulas = label.split(',').slice(1).map(x => x.trim().toLowerCase());
    if (names === lower || formulas.includes(lower)) {
      return { mol: moleculeFromSMILES(smiles), label: label, source: 'database' };
    }
  }
  // 2. try as SMILES (contains =, #, (, ), digits attached to atoms, ring digits, or pure element symbols)
  const looksLikeSMILES = /[=#@]/.test(q) || /^c|^n|^o|^s|^p|^f|^cl|^br/i.test(q);
  if (looksLikeSMILES) {
    try {
      return { mol: moleculeFromSMILES(q), label: q, source: 'SMILES' };
    } catch (e) {
      // fall through to formula attempt
      var smileErr = e;
    }
  }
  // 3. try as molecular formula: parse, then generate a connectivity guess
  // For formulas we can only guess; use a heuristic builder for simple known ones,
  // otherwise report the formula composition.
  try {
    const counts = parseFormula(q);
    return { formulaCounts: counts, label: formulaString(counts), source: 'formula' };
  } catch (e) {
    throw smileErr || new Error('Could not interpret input: ' + q);
  }
}

// Formula -> approximate structure for common simple cases via template table.
// Many formulas map to a unique small molecule; for others the app reports composition.
const FORMULA_TEMPLATES = {
  H2O: 'O', CH4: 'C', NH3: 'N', HCl: 'Cl', CO2: 'O=C=O', O2: 'O=O', N2: 'N#N',
  H2O2: 'OO', CH3OH: 'CO', C2H6: 'CC', C2H4: 'C=C', C2H2: 'C#C', CH2O: 'C=O',
  C2H6O: 'CCO', CH4O: 'CO', C3H8: 'CCC', C4H10: 'CCCC', C6H6: 'c1ccccc1',
  C6H12O6: 'OCC1OC(O)C(O)C(O)C1O', NaCl: '[Na+].[Cl-]', NH4: '[NH4+]'
};

function moleculeFromFormula(f) {
  const counts = parseFormula(f);
  const key = formulaString(counts);
  if (FORMULA_TEMPLATES[key]) return moleculeFromSMILES(FORMULA_TEMPLATES[key]);
  return null;
}
