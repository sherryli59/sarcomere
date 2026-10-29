// Execute the saved viewer against a minimal DOM/canvas to check every case
// and realization. This checks code/data wiring, not browser rendering.
'use strict';
const fs = require('node:fs');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const html = fs.readFileSync(process.argv[2], 'utf8');
const payload = html.match(/<script id="data" type="application\/json">([\s\S]*?)<\/script>/)[1];
const code = html.match(/<\/script><script>([\s\S]*?)<\/script>/)[1];
const elements = new Map();
let calls = 0;
const drawing = new Proxy({}, {
  get(target, key) {
    if (key in target) return target[key];
    return (...args) => {
      for (const value of args) {
        if (typeof value === 'number') assert(Number.isFinite(value), `Nonfinite canvas ${key}`);
      }
      calls++;
    };
  },
});
function element(id) {
  if (!elements.has(id)) elements.set(id, {
    value: '0', checked: true, textContent: '', innerHTML: '',
    getContext: () => drawing,
    getBoundingClientRect: () => ({width: 1200, height: 460}),
    setPointerCapture() {},
  });
  return elements.get(id);
}
element('data').textContent = payload;
const sandbox = {
  document: {getElementById: element},
  window: {devicePixelRatio: 1, addEventListener() {}},
  requestAnimationFrame() {},
};
vm.createContext(sandbox);
vm.runInContext(code, sandbox);
const data = sandbox.window.standaloneContact.D;
let frames = 0, realizations = 0;
for (let i = 0; i < data.cases.length; i++) {
  const c = data.cases[i];
  element('case').value = String(i);
  vm.runInContext('experiment()', sandbox);
  for (let rep = 0; rep < c.center[0].length; rep++) {
    realizations++;
    element('rep').value = String(rep);
    for (const frame of [0, Math.floor(c.time.length/2), c.time.length-1]) {
      element('time').value = String(frame);
      sandbox.window.standaloneContact.draw();
      assert(element('readout').innerHTML.includes(`realization ${rep} /`));
      assert(!/NaN|undefined|Infinity/.test(element('readout').innerHTML));
      const live = c.alive[frame][rep];
      for (let j = 0; j < live.length; j++) {
        if (!live[j]) {
          assert(c.fi[frame][rep][j].every(x => x === 0));
          assert(c.fj[frame][rep][j].every(x => x === 0));
        }
      }
      frames++;
    }
  }
  for (const id of ['arrows','anchors','dead']) {
    element(id).checked = false;
    element(id).onchange();
    element(id).checked = true;
  }
  element('view').onclick();
  element('reset').onclick();
  element('play').onclick();
  element('play').onclick();
}
assert.equal(data.cases.length, 6);
assert.equal((html.match(/data:image\/png;base64,/g) || []).length, 4);
assert(!html.includes('__DATA__') && !html.includes('__FIGURES__'));
console.log(JSON.stringify({passed:true, cases:data.cases.length,realizations,frames,
  canvas_calls:calls,scope:'DOM/canvas wiring smoke test; actual browser rendering not tested'}));
