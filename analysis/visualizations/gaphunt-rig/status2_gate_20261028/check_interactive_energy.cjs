// Cross-check browser math against independent Python outputs and exercise UI wiring.
// Plotly rendering requires an actual browser; this smoke check mocks only its API.
const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
const html=fs.readFileSync(process.argv[2],'utf8');
const elements=new Map();
function element(id){if(!elements.has(id))elements.set(id,{id,checked:true,value:'',style:{},dataset:{},textContent:'',innerHTML:'',querySelectorAll:()=>[],on:()=>{}});return elements.get(id);}
const checks=['body','boundary','approach','slice','examples','forces','endpoint','reaction','shell-0','shell-1','shell-2','shell-3'];
checks.forEach(element);['approach','slice','endpoint','reaction'].forEach(id=>element(id).checked=false);
['px','py','pz'].forEach((id,i)=>element(id).value=[750,80,0][i]);
element('mode').value='0';element('opacity').value='.16';element('scale').value='4';
element('energy-data').textContent=html.match(/<script id="energy-data"[^>]*>([^<]+)<\/script>/)[1];
const views=['whole','tip','cross'].map(view=>({dataset:{view}}));
let traces=[],layout;
const context={console,Uint8Array,atob,Blob,Response,DecompressionStream,setTimeout,
 document:{getElementById:element,body:{dataset:{}},querySelectorAll:selector=>selector.includes('checkbox')?checks.map(element):selector.includes('data-view')?views:[]},
 Plotly:{newPlot:(id,t,l)=>{traces=t;layout=l;return Promise.resolve();},react:(id,t,l)=>{traces=t;layout=l;return Promise.resolve();}},
 window:{addEventListener:()=>{}}};
vm.createContext(context);
const scripts=[...html.matchAll(/<script(?:\s[^>]*)?>([\s\S]*?)<\/script>/g)];
vm.runInContext(scripts.at(-1)[1],context);
(async()=>{
 for(let n=0;n<100&&context.document.body.dataset.ready!=='true';n++)await new Promise(r=>setTimeout(r,30));
 assert.equal(context.document.body.dataset.ready,'true');
 const {evaluate,DATA}=context.window.energyViewer;
 let maxEnergyError=0,maxForceError=0;
 for(let m=0;m<DATA.modes.length;m++)for(let i=0;i<DATA.examples.length;i++){
  const v=evaluate(DATA.examples[i].position,m),reference=DATA.modes[m];
  maxEnergyError=Math.max(maxEnergyError,Math.abs(v.energy-reference.referenceEnergy[i]));
  v.f.forEach((x,j)=>maxForceError=Math.max(maxForceError,Math.abs(x-reference.referenceForce[i][j])));
  assert.ok(Math.hypot(...v.f)<=30+1e-9);
  v.f.forEach((x,j)=>assert.ok(Math.abs(x-v.fi[j]-v.fj[j])<1e-10));
 }
 assert.ok(maxEnergyError<1e-5);assert.ok(maxForceError<1e-9);
 const full=traces.length;element('shell-0').checked=false;element('shell-0').onchange();assert.equal(traces.length,full-1);
 element('forces').checked=false;element('forces').onchange();const noForces=traces.length;
 element('forces').checked=true;element('forces').onchange();assert.ok(traces.length>noForces);
 for(let m=0;m<4;m++){element('mode').value=String(m);element('mode').onchange();assert.ok(element('mode-note').textContent.length>0);}
 element('slice').checked=true;element('slice').onchange();assert.ok(traces.some(t=>t.type==='surface'));
 element('px').value='700';element('py').value='140';element('pz').value='0';element('py').onchange();
 assert.ok(element('probe-readout').innerHTML.includes('Custom probe'));
 const previous=traces.length;element('endpoint').checked=true;element('endpoint').onchange();assert.ok(traces.length>previous);
 views[1].onclick();assert.equal(layout.scene.xaxis.range[0],580);
 element('scale').value='8';element('scale').oninput();assert.equal(element('scale-value').textContent,8);
 fs.writeFileSync(process.argv[3],JSON.stringify({passed:true,browser_rendering_checked:false,max_python_js_energy_error_pNnm:maxEnergyError,max_python_js_force_error_pN:maxForceError,
  checked:'4 modes; 8 examples; force bound; endpoint reactions; contour, force, slice and endpoint toggles; custom probe; view and scale controls'},null,2));
 console.log('Browser math / control smoke checks passed.',maxEnergyError,maxForceError);
})().catch(e=>{console.error(e);process.exitCode=1;});
