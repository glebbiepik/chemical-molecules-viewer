// viewer.js — three.js rendering (ball-and-stick / space-fill / wireframe)

function createViewer(container) {
  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x0b0e14);

  const camera = new THREE.PerspectiveCamera(45, 1, 0.1, 2000);
  const renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  container.appendChild(renderer.domElement);

  // lighting
  scene.add(new THREE.AmbientLight(0xffffff, 0.45));
  const key = new THREE.DirectionalLight(0xffffff, 1.1);
  key.position.set(5, 8, 10);
  scene.add(key);
  const fill = new THREE.DirectionalLight(0x8899ff, 0.35);
  fill.position.set(-6, -4, -8);
  scene.add(fill);

  const moleculeGroup = new THREE.Group();
  scene.add(moleculeGroup);

  // simple orbit controls (no external OrbitControls file needed)
  const controls = {
    target: new THREE.Vector3(),
    distance: 20,
    theta: 0.6, phi: 1.1,
    dragging: false, lastX: 0, lastY: 0, panning: false
  };

  function updateCamera() {
    const { target, distance, theta, phi } = controls;
    const sp = Math.max(0.05, Math.min(Math.PI - 0.05, phi));
    camera.position.set(
      target.x + distance * Math.sin(sp) * Math.sin(theta),
      target.y + distance * Math.cos(sp),
      target.z + distance * Math.sin(sp) * Math.cos(theta)
    );
    camera.lookAt(target);
  }

  const dom = renderer.domElement;
  dom.style.touchAction = 'none';
  dom.addEventListener('pointerdown', e => {
    controls.dragging = true;
    controls.panning = e.button === 2 || e.shiftKey;
    controls.lastX = e.clientX; controls.lastY = e.clientY;
    dom.setPointerCapture(e.pointerId);
  });
  dom.addEventListener('pointermove', e => {
    if (!controls.dragging) return;
    const dx = e.clientX - controls.lastX, dy = e.clientY - controls.lastY;
    controls.lastX = e.clientX; controls.lastY = e.clientY;
    if (controls.panning) {
      const panScale = controls.distance * 0.0016;
      const right = new THREE.Vector3().setFromMatrixColumn(camera.matrix, 0);
      const up = new THREE.Vector3().setFromMatrixColumn(camera.matrix, 1);
      controls.target.addScaledVector(right, -dx * panScale);
      controls.target.addScaledVector(up, dy * panScale);
    } else {
      controls.theta -= dx * 0.006;
      controls.phi -= dy * 0.006;
      controls.phi = Math.max(0.05, Math.min(Math.PI - 0.05, controls.phi));
    }
    updateCamera();
  });
  dom.addEventListener('pointerup', e => { controls.dragging = false; dom.releasePointerCapture(e.pointerId); });
  dom.addEventListener('contextmenu', e => e.preventDefault());
  dom.addEventListener('wheel', e => {
    e.preventDefault();
    controls.distance = Math.max(2, Math.min(300, controls.distance * (1 + Math.sign(e.deltaY) * 0.1)));
    updateCamera();
  }, { passive: false });

  function resize() {
    const w = container.clientWidth, h = container.clientHeight;
    renderer.setSize(w, h);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
  }
  window.addEventListener('resize', resize);

  let spin = true;
  function animate() {
    requestAnimationFrame(animate);
    if (spin && !controls.dragging) { controls.theta += 0.0035; updateCamera(); }
    renderer.render(scene, camera);
  }

  // ---------- molecule building ----------
  const SPHERE_SEG = 32;
  function build(mol, opts) {
    const style = opts.style || 'ball-stick';
    const scale = opts.scale || 1.0;
    while (moleculeGroup.children.length) {
      const c = moleculeGroup.children.pop();
      if (c.geometry) c.geometry.dispose();
      if (c.material) c.material.dispose();
    }
    if (!mol || !mol.pos || !mol.pos.length) return;

    // bounds & centering
    let min = [1e9, 1e9, 1e9], max = [-1e9, -1e9, -1e9];
    for (const p of mol.pos) for (let k = 0; k < 3; k++) { min[k] = Math.min(min[k], p[k ? ['x','y','z'][k] : 'x']); }
    for (const p of mol.pos) {
      min[0] = Math.min(min[0], p.x); max[0] = Math.max(max[0], p.x);
      min[1] = Math.min(min[1], p.y); max[1] = Math.max(max[1], p.y);
      min[2] = Math.min(min[2], p.z); max[2] = Math.max(max[2], p.z);
    }
    const center = new THREE.Vector3((min[0]+max[0])/2, (min[1]+max[1])/2, (min[2]+max[2])/2);
    const size = Math.max(max[0]-min[0], max[1]-min[1], max[2]-min[2]) || 4;

    const radiusScale = style === 'space-fill' ? 2.6 : 0.55;
    const sphereGeoms = new Map();
    function sphereGeom(r) {
      const key = Math.round(r * 100);
      if (!sphereGeoms.has(key)) sphereGeoms.set(key, new THREE.SphereGeometry(r, SPHERE_SEG, Math.max(12, SPHERE_SEG / 2)));
      return sphereGeoms.get(key);
    }
    const cylGeom = new THREE.CylinderGeometry(1, 1, 1, 24, 1, true);

    for (let i = 0; i < mol.atoms.length; i++) {
      const a = mol.atoms[i];
      const el = ELEMENTS[a.el] || { color: 0xff00ff, rcov: 1.0 };
      const p = new THREE.Vector3(mol.pos[i].x, mol.pos[i].y, mol.pos[i].z).sub(center).multiplyScalar(scale);
      let r;
      if (style === 'space-fill') r = Math.max(1.1, el.rcov * 2.2);
      else r = Math.max(0.28, el.rcov * (a.el === 'H' ? 0.85 : 0.62));
      const mesh = new THREE.Mesh(sphereGeom(r * scale), new THREE.MeshPhongMaterial({ color: el.color, shininess: 60, specular: 0x333333 }));
      mesh.position.copy(p);
      moleculeGroup.add(mesh);
    }

    if (style !== 'space-fill') {
      for (const b of mol.bonds) {
        const pa = new THREE.Vector3(mol.pos[b.i].x, mol.pos[b.i].y, mol.pos[b.i].z).sub(center).multiplyScalar(scale);
        const pb = new THREE.Vector3(mol.pos[b.j].x, mol.pos[b.j].y, mol.pos[b.j].z).sub(center).multiplyScalar(scale);
        const elA = ELEMENTS[mol.atoms[b.i].el], elB = ELEMENTS[mol.atoms[b.j].el];
        const cA = new THREE.Color(elA ? elA.color : 0xffffff);
        const cB = new THREE.Color(elB ? elB.color : 0xffffff);
        const radius = style === 'wireframe' ? 0.035 * scale : 0.09 * scale;
        const orders = (mol.atoms[b.i].aromatic && mol.atoms[b.j].aromatic) ? [1, 2] : b.order === 3 ? [1, 1, 1] : b.order === 2 ? [1, 1] : [1];
        // spread multiple bond cylinders perpendicular to bond dir
        const dir = new THREE.Vector3().subVectors(pb, pa);
        const len = dir.length();
        dir.normalize();
        let perp = new THREE.Vector3(0, 1, 0).cross(dir);
        if (perp.lengthSq() < 0.01) perp = new THREE.Vector3(1, 0, 0).cross(dir);
        perp.normalize();
        const spread = 0.13 * scale;
        orders.forEach((_, k) => {
          const offset = orders.length === 1 ? new THREE.Vector3() : perp.clone().multiplyScalar((k - (orders.length - 1) / 2) * spread);
          const mid = new THREE.Vector3().addVectors(pa, pb).multiplyScalar(0.5).add(offset);
          const a1 = pa.clone().add(offset), b1 = pb.clone().add(offset);
          // two half-bonds colored by atom
          const segs = [[a1, mid, cA], [mid, b1, cB]];
          for (const [s0, s1, col] of segs) {
            const segLen = s0.distanceTo(s1);
            if (segLen < 1e-4) continue;
            const mesh = new THREE.Mesh(cylGeom, new THREE.MeshPhongMaterial({ color: col, shininess: 40 }));
            mesh.scale.set(radius, segLen, radius);
            mesh.position.copy(s0).add(s1).multiplyScalar(0.5);
            mesh.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), new THREE.Vector3().subVectors(s1, s0).normalize());
            moleculeGroup.add(mesh);
          }
        });
      }
    }

    // fit camera
    controls.target.set(0, 0, 0);
    controls.distance = size * scale * 1.9 + 4;
    updateCamera();
    resize();
  }

  animate();
  resize();
  updateCamera();

  return {
    build,
    setSpin(v) { spin = v; },
    resetView() { controls.theta = 0.6; controls.phi = 1.1; controls.target.set(0, 0, 0); updateCamera(); },
    resize
  };
}
