// app.js — UI glue

const inputEl = document.getElementById('input');
const goBtn = document.getElementById('goBtn');
const errorEl = document.getElementById('error');
const infoEl = document.getElementById('info');
const styleEl = document.getElementById('style');
const spinEl = document.getElementById('spin');
const listEl = document.getElementById('molList');

const viewer = createViewer(document.getElementById('viewport'));
let current = null;   // { mol, label, source }

function showError(msg) { errorEl.textContent = msg || ''; }

function showInfo(r) {
  const mol = r.mol;
  const counts = {};
  for (const a of mol.atoms) counts[a.el] = (counts[a.el] || 0) + 1;
  const comp = Object.keys(counts).map(k => k + ':' + counts[k]).join(' · ');
  infoEl.style.display = 'block';
  infoEl.innerHTML =
    '<b>' + escapeHtml(r.label) + '</b><br>' +
    'Formula: ' + escapeHtml(mol.formula) + '<br>' +
    'Molar mass: ' + mol.molarMass.toFixed(2) + ' g/mol<br>' +
    'Atoms: ' + mol.atoms.length + ' · Bonds: ' + mol.bonds.length + '<br>' +
    '<span style="color:#7a8199">' + escapeHtml(comp) + '</span><br>' +
    '<span style="color:#5b9dff;font-size:0.72rem">via ' + escapeHtml(r.source) + '</span>';
}

function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
}

function build(query) {
  showError('');
  try {
    const r = resolveInput(query);
    if (!r.mol) {
      // formula without a template: report composition and suggest
      const comp = Object.keys(r.formulaCounts).map(k => k + ':' + r.formulaCounts[k]).join(' · ');
      showError('Structure for "' + r.label + '" is ambiguous. Composition: ' + comp + '. Try a SMILES string or a name from the list.');
      return;
    }
    current = r;
    viewer.build(r.mol, { style: styleEl.value, scale: 1 });
    showInfo(r);
  } catch (e) {
    showError(e.message);
  }
}

goBtn.addEventListener('click', () => build(inputEl.value));
inputEl.addEventListener('keydown', e => { if (e.key === 'Enter') build(inputEl.value); });
styleEl.addEventListener('change', () => { if (current) viewer.build(current.mol, { style: styleEl.value, scale: 1 }); });
spinEl.addEventListener('change', () => viewer.setSpin(spinEl.checked));

// quick picks
const picks = ['water', 'methane', 'ethanol', 'benzene', 'caffeine', 'glucose', 'aspirin', 'testosterone', 'CO2', 'nh3'];
for (const name of picks) {
  const chip = document.createElement('div');
  chip.className = 'chip';
  const dbEntry = DB[keyOf(name)];
  chip.textContent = dbEntry ? dbEntry[0].split(',')[0] : name;
  chip.addEventListener('click', () => { inputEl.value = name; build(name); });
  listEl.appendChild(chip);
}
function keyOf(name) {
  for (const k in DB) if (k === name.toLowerCase().replace(/[\s-]/g, '_')) return k;
  return name.toLowerCase().replace(/[\s-]/g, '_');
}

// initial molecule
build('water');
