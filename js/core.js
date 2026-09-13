// core.js — chemistry engine: element data, formula/SMILES parsing, 3D geometry
// No dependencies. Coordinates in Angstroms.

const ELEMENTS = {
  H:  { z:1,  m:1.008,  rcov:0.31, color:0xf8f8f8, name:'Hydrogen'  },
  B:  { z:5,  m:10.81,  rcov:0.84, color:0xffb5b5, name:'Boron'     },
  C:  { z:6,  m:12.011, rcov:0.76, color:0x909090, name:'Carbon'    },
  N:  { z:7,  m:14.007, rcov:0.71, color:0x3050f8, name:'Nitrogen'  },
  O:  { z:8,  m:15.999, rcov:0.66, color:0xff4d4d, name:'Oxygen'    },
  F:  { z:9,  m:18.998, rcov:0.57, color:0x90e050, name:'Fluorine'  },
  Si: { z:14, m:28.085, rcov:1.11, color:0xf0c8a0, name:'Silicon'   },
  P:  { z:15, m:30.974, rcov:1.07, color:0xff8000, name:'Phosphorus'},
  S:  { z:16, m:32.06,  rcov:1.05, color:0xffd52a, name:'Sulfur'    },
  Cl: { z:17, m:35.45,  rcov:1.02, color:0x1ff01f, name:'Chlorine'  },
  Br: { z:35, m:79.904, rcov:1.20, color:0xa62929, name:'Bromine'   },
  I:  { z:53, m:126.9,  rcov:1.39, color:0x940094, name:'Iodine'    }
};

// Ideal bond lengths (Angstrom) by atom pair — geometric defaults
const BOND_LEN = {
  'C-C': 1.54, 'C=C': 1.34, 'C#C': 1.20, 'C~C': 1.40,
  'C-N': 1.47, 'C=N': 1.28, 'C#N': 1.16, 'C~N': 1.34,
  'C-O': 1.43, 'C=O': 1.23, 'C~O': 1.30,
  'C-H': 1.09, 'N-H': 1.01, 'O-H': 0.96,
  'C-F': 1.35, 'C-Cl': 1.77, 'C-Br': 1.94, 'C-I': 2.14,
  'C-S': 1.82, 'C-P': 1.84,
  'N-N': 1.45, 'N-O': 1.40, 'N=O': 1.21,
  'O-O': 1.48, 'S-O': 1.58, 'P-O': 1.61, 'S-H': 1.34, 'P=O': 1.48,
  'N~N': 1.35, 'N~O': 1.36
};

function bondLength(a, b, order, aromatic) {
  const key = a + (order === 2 ? '=' : order === 3 ? '#' : aromatic ? '~' : '-') + b;
  if (BOND_LEN[key]) return BOND_LEN[key];
  const rev = b + (order === 2 ? '=' : order === 3 ? '#' : aromatic ? '~' : '-') + a;
  if (BOND_LEN[rev]) return BOND_LEN[rev];
  const ea = ELEMENTS[a], eb = ELEMENTS[b];
  if (!ea || !eb) return 1.5;
  let base = (ea.rcov + eb.rcov) / 0.728; // convert cov radii sum to ~bond length
  if (order === 2) base *= 0.91;
  if (order === 3) base *= 0.85;
  return base;
}

class Molecule {
  constructor() { this.atoms = []; this.bonds = []; }
  addAtom(el, aromatic = false) {
    const a = { el, aromatic, x: 0, y: 0, z: 0 };
    this.atoms.push(a);
    return this.atoms.length - 1;
  }
  addBond(i, j, order = 1) {
    this.bonds.push({ i, j, order });
  }
  neighbors(i) {
    const out = [];
    for (const b of this.bonds) {
      if (b.i === i) out.push({ j: b.j, order: b.order, bond: b });
      else if (b.j === i) out.push({ j: b.i, order: b.order, bond: b });
    }
    return out;
  }
  bondBetween(i, j) {
    return this.bonds.find(b => (b.i === i && b.j === j) || (b.i === j && b.j === i));
  }
  bondOrderSum(i) {
    return this.neighbors(i).reduce((s, n) => s + n.order, 0);
  }
}

// ---------- Formula parser ----------
function parseFormula(f) {
  const s = f.replace(/\s+/g, '');
  let i = 0;
  function group() {
    const counts = {};
    while (i < s.length) {
      const c = s[i];
      if (c === '(') {
        i++;
        const sub = group();
        if (s[i] !== ')') throw new Error('Mismatched parentheses');
        i++;
        let n = number();
        for (const k in sub) counts[k] = (counts[k] || 0) + sub[k] * n;
      } else if (c === ')') {
        return counts;
      } else if (/[A-Z]/.test(c)) {
        let sym = c;
        i++;
        if (i < s.length && /[a-z]/.test(s[i])) sym += s[i++];
        const n = number();
        counts[sym] = (counts[sym] || 0) + n;
      } else {
        throw new Error('Unexpected character "' + c + '"');
      }
    }
    return counts;
  }
  function number() {
    let num = '';
    while (i < s.length && /[0-9]/.test(s[i])) num += s[i++];
    return num ? parseInt(num, 10) : 1;
  }
  const counts = group();
  if (!Object.keys(counts).length) throw new Error('Empty formula');
  return counts;
}

function formulaString(counts) {
  const keys = Object.keys(counts);
  const order = [];
  if (keys.includes('C')) { order.push('C'); if (keys.includes('H')) order.push('H'); }
  keys.sort();
  for (const k of keys) if (!order.includes(k)) order.push(k);
  return order.map(k => k + (counts[k] > 1 ? counts[k] : '')).join('');
}

// ---------- Molecule database ----------
const DB = {
  water:        ['water, H2O',                        'O'],
  methane:      ['methane, CH4',                      'C'],
  ammonia:      ['ammonia, NH3',                      'N'],
  hcl:          ['hydrogen chloride, HCl',            'Cl'],
  co2:          ['carbon dioxide, CO2',               'O=C=O'],
  hydrogen_peroxide: ['hydrogen peroxide, H2O2',      'OO'],
  methanol:     ['methanol, CH3OH',                   'CO'],
  ethanol:      ['ethanol, C2H5OH',                   'CCO'],
  dimethyl_ether: ['dimethyl ether, C2H6O',           'COC'],
  formaldehyde: ['formaldehyde, CH2O',                'C=O'],
  formic_acid:  ['formic acid, HCOOH',                'OC=O'],
  acetic_acid:  ['acetic acid, CH3COOH',              'CC(=O)O'],
  acetone:      ['acetone, C3H6O',                    'CC(=O)C'],
  lactic_acid:  ['lactic acid, C3H6O3',               'CC(O)C(=O)O'],
  urea:         ['urea, CH4N2O',                      'NC(=O)N'],
  glycine:      ['glycine, C2H5NO2',                  'NCC(=O)O'],
  alanine:      ['alanine, C3H7NO2',                  'CC(N)C(=O)O'],
  ethane:       ['ethane, C2H6',                      'CC'],
  ethene:       ['ethene (ethylene), C2H4',           'C=C'],
  ethyne:       ['ethyne (acetylene), C2H2',          'C#C'],
  propane:      ['propane, C3H8',                     'CCC'],
  butane:       ['butane, C4H10',                     'CCCC'],
  isobutane:    ['isobutane, C4H10',                  'C(C)(C)C'],
  cyclohexane:  ['cyclohexane, C6H12',                'C1CCCCC1'],
  benzene:      ['benzene, C6H6',                     'c1ccccc1'],
  toluene:      ['toluene, C7H8',                     'Cc1ccccc1'],
  phenol:       ['phenol, C6H5OH',                    'Oc1ccccc1'],
  aniline:      ['aniline, C6H5NH2',                  'Nc1ccccc1'],
  styrene:      ['styrene, C8H8',                     'C=Cc1ccccc1'],
  naphthalene:  ['naphthalene, C10H8',                'c1ccc2ccccc2c1'],
  glucose:      ['glucose (ring), C6H12O6',           'OCC1OC(O)C(O)C(O)C1O'],
  sucrose:      ['sucrose, C12H22O11',                'OCC1OC(OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O'],
  aspirin:      ['aspirin, C9H8O4',                   'CC(=O)Oc1ccccc1C(=O)O'],
  caffeine:     ['caffeine, C8H10N4O2',               'Cn1c(=O)c2c(ncn2C)n(C)c1=O'],
  limonene:     ['limonene, C10H16',                  'CC(=C)C1CCC(C)=CCC1'],
  testosterone: ['testosterone, C19H28O2',            'CC12CCC3C(CCC(=O)C3=C1)CCC1C2CCC(O)C1'],
  cholesterol:  ['cholesterol, C27H46O',              'CC(C)CCCC(C)C1CCC2(C)C1CCC1C2CC=C2CC(O)CCC12C'],
  penicillin:   ['penicillin G, C16H18N2O4S',         'CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O'],
  atp:          ['ATP, C10H16N5O13P3',                'Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O']
};

// ---------- SMILES parser (organic subset) ----------
function parseSMILES(smiles) {
  const s = smiles.trim();
  if (!s) throw new Error('Empty SMILES');
  const mol = new Molecule();
  let i = 0;
  const branchStack = [];   // {idx, order}
  let prev = null;
  let pendingOrder = 1;
  const rings = {};         // digit -> {atom, order}
  const aromaticSet = new Set();

  function readBondSymbol() {
    let o = pendingOrder;
    pendingOrder = 1;
    while (i < s.length) {
      if (s[i] === '=') { o = 2; i++; }
      else if (s[i] === '#') { o = 3; i++; }
      else if (s[i] === '/' || s[i] === '\\') { i++; if (o === 1) o = 1; }
      else break;
    }
    return o;
  }

  function readRingDigit() {
    let num = '';
    while (i < s.length && /[0-9]/.test(s[i])) num += s[i++];
    return num ? parseInt(num, 10) : null;
  }

  while (i < s.length) {
    let c = s[i];
    if (c === '(') { i++; branchStack.push({ idx: prev, order: pendingOrder }); pendingOrder = 1; continue; }
    if (c === ')') {
      i++;
      if (!branchStack.length) throw new Error('Unmatched ")"');
      const p = branchStack.pop();
      prev = p.idx; pendingOrder = p.order;
      continue;
    }
    if (c === '.') { i++; prev = null; pendingOrder = 1; continue; }

    let order = 1;
    if (prev !== null && (c === '=' || c === '#' || c === '/' || c === '\\')) {
      order = readBondSymbol();
      c = s[i]; // re-read after consuming bond symbol(s)
    }

    // ring closure?
    if (/[0-9]/.test(c)) {
      // ring digits must not swallow element symbols that follow (e.g. C12... is C, ring 1, ring 2 handled greedily only if valid)
      if (prev === null) throw new Error('Ring digit without preceding atom');
      // A digit run is a two-digit ring number only if that ring is already
      // open; otherwise each digit opens/closes its own ring (SMILES "C12").
      let numStr = '';
      while (i < s.length && /[0-9]/.test(s[i])) { numStr += s[i]; i++; }
      const nums = [];
      if (rings[parseInt(numStr, 10)] !== undefined) {
        nums.push(parseInt(numStr, 10));
      } else {
        for (const ch of numStr) nums.push(parseInt(ch, 10));
      }
      for (const num of nums) {
        if (rings[num] !== undefined) {
          mol.addBond(rings[num].atom, prev, rings[num].order !== 1 ? rings[num].order : (order !== 1 ? order : 1));
          delete rings[num];
        } else {
          rings[num] = { atom: prev, order };
        }
      }
      continue;
    }

    let el, aromatic = false;
    if (c === '[') {
      i++;
      let tok = '';
      while (i < s.length && s[i] !== ']') { tok += s[i]; i++; }
      i++;
      tok = tok.replace(/@/g, '');
      const m = tok.match(/^(se|as|[A-Z][a-z]?)/);
      if (!m) throw new Error('Cannot read element from "' + tok + '"');
      el = m[1];
      if (el[0] === el[0].toLowerCase()) { aromatic = true; el = el.toUpperCase(); }
      if (!ELEMENTS[el]) throw new Error('Unknown element "' + el + '"');
    } else if (/[A-Z]/.test(c)) {
      let sym = c; i++;
      if (i < s.length && /[a-z]/.test(s[i])) {
        const two = c + s[i];
        if (ELEMENTS[two]) { sym = two; i++; }
      }
      el = sym;
    } else if (/[bcnosp]/.test(c)) {
      el = c.toUpperCase(); aromatic = true; i++;
    } else {
      throw new Error('Unexpected character "' + c + '" at position ' + i);
    }

    const idx = mol.addAtom(el, aromatic);
    if (aromatic) aromaticSet.add(idx);
    if (prev !== null) {
      let o = order;
      if (aromatic && mol.atoms[prev].aromatic && o === 1) o = 1; // keep 1; kekulize later
      mol.addBond(prev, idx, o);
    }
    prev = idx;
  }
  if (Object.keys(rings).length) throw new Error('Unclosed ring digit(s): ' + Object.keys(rings).join(','));
  kekulize(mol, aromaticSet);
  inferHydrogens(mol);
  return mol;
}

// Convert aromatic (lowercase) ring bonds to alternating 1/2 orders.
// Bonds are assigned in creation order, which follows the ring path in SMILES,
// so alternation yields a valid Kekule structure for typical rings.
function kekulize(mol, aromaticSet) {
  if (!aromaticSet.size) return;
  let toggle = 1;
  for (const b of mol.bonds) {
    if (aromaticSet.has(b.i) && aromaticSet.has(b.j)) {
      b.order = toggle;
      toggle = toggle === 1 ? 2 : 1;
    }
  }
}

// Attach implicit hydrogens based on standard valences.
function inferHydrogens(mol) {
  const VALENCE = { C: 4, N: 3, O: 2, S: 2, P: 3, F: 1, Cl: 1, Br: 1, I: 1, B: 3, Si: 4 };
  const toAdd = [];
  for (let idx = 0; idx < mol.atoms.length; idx++) {
    const a = mol.atoms[idx];
    const v = VALENCE[a.el];
    if (!v) continue;
    let bos = mol.bondOrderSum(idx);
    let need = v - bos;
    if (need < 0) need = 0; // charged / hypervalent: skip
    for (let k = 0; k < need; k++) toAdd.push(idx);
  }
  for (const idx of toAdd) {
    const h = mol.addAtom('H');
    mol.addBond(idx, h, 1);
  }
}
