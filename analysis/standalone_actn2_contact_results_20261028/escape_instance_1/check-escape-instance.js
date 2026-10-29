const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
const path=process.argv[2],html=fs.readFileSync(path,'utf8');
const payload=html.match(/<script type="application\/json" id="myosin-escape-one-data">([\s\S]*?)<\/script>/)[1];
const code=html.match(/<script>\s*([\s\S]*?)<\/script>/)[1];
let width=736;
class Element{
  constructor(tag='div'){this.tag=tag;this.attrs={};this.children=[];this.listeners={};this.textContent='';this.value='0';}
  setAttribute(k,v){this.attrs[k]=v;}
  getAttribute(k){return this.attrs[k];}
  appendChild(c){this.children.push(c);return c;}
  replaceChildren(){this.children=[];}
  addEventListener(k,v){this.listeners[k]=v;}
  getBoundingClientRect(){return {width};}
}
const root=new Element(),svg=new Element('svg'),nodes=new Map();
for(const key of ['[data-time-slider]','[data-play]','[data-status]','[data-time]','[data-announcement]'])nodes.set(key,new Element());
nodes.set('.escape-drawing',svg);
root.querySelector=s=>{assert(nodes.has(s),s);return nodes.get(s);};
const buttons=['before','after','escape'].map(x=>{const e=new Element('button');e.setAttribute('data-step',x);return e;});
root.querySelectorAll=s=>{assert.equal(s,'[data-step]');return buttons;};
const sandbox={document:{getElementById:id=>id==='myosin-escape-one'?root:{textContent:payload},createElementNS:(ns,tag)=>new Element(tag)},window:{addEventListener(){}},ResizeObserver:class{constructor(f){this.f=f;}observe(){this.f();}},requestAnimationFrame:()=>1,cancelAnimationFrame(){},performance:{now:()=>0}};
vm.createContext(sandbox);vm.runInContext(code,sandbox);
const viewer=sandbox.window.myosinEscapeInstance;
assert.equal(viewer.stages.before,11);assert.equal(viewer.stages.after,12);
const snapshots=[];
for(width of [320,360,736,1024]){
  for(const i of [0,11,12,20,40,80,200]){
    viewer.select(i);
    assert(!/NaN|undefined|Infinity/.test(nodes.get('[data-status]').textContent));
    const dims=svg.attrs.viewBox.split(' ').map(Number);
    for(const child of svg.children){
      for(const [k,v] of Object.entries(child.attrs))assert(!/NaN|undefined|Infinity/.test(String(v)),`${child.tag} ${k}`);
      if(child.tag==='text'){
        const x=Number(child.attrs.x),y=Number(child.attrs.y);
        assert(x>=0&&x<=dims[2]&&y>=0&&y<=dims[3],`Text outside viewport ${child.textContent}`);
      }
    }
    if(width===736)snapshots.push({frame:i,status:nodes.get('[data-status]').textContent});
  }
}
for(const b of buttons)b.listeners.click();
assert.equal(nodes.get('[data-time-slider]').value,'200');
const data=viewer.data;
assert(data.alive[11].every(Boolean));assert.deepEqual(Array.from(data.alive[12]),[0,1,1,1]);
assert.deepEqual(Array.from(data.alive[200]),[0,1,1,1]);
assert(Math.hypot(...data.center[11].slice(1))<1e-6);
assert(Math.hypot(...data.center[200].slice(1))>102);
assert(data.center[200][0]>1067);
assert(html.length<1000000);
console.log(JSON.stringify({passed:true,viewport_widths:[320,360,736,1024],snapshots,scope:'Saved data and SVG/DOM wiring; browser rendering not tested'},null,2));
