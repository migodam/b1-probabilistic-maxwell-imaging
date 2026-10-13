// Code-level verification only. This does not launch a browser or render the page.
const fs = require('fs');
const path = require('path');
const vm = require('vm');
const assert = require('assert');
const started = process.hrtime.bigint();
const cpu0 = process.cpuUsage();
const root = path.resolve(__dirname, '../../..');
const html = fs.readFileSync(path.join(root, 'dashboard.html'), 'utf8');
assert(!/__BLOB__|__ORDER__/.test(html), 'unexpanded data placeholders');
assert(!/<(?:script|link|iframe)[^>]*(?:src|href)\s*=/i.test(html), 'external script/style/frame reference');
assert(!/<img[^>]*src\s*=\s*["']https?:/i.test(html), 'remote image reference');
assert(/@media\(max-width:800px\)/.test(html), 'responsive rules absent');
const script = html.match(/<script>([\s\S]*?)<\/script>/)[1];
class Element {
  constructor(id='') { this.id=id; this.dataset={}; this.innerHTML=''; this.textContent=''; this.src=''; this.children=[]; this.classes=new Set(); this.events={}; }
  get classList() { return {toggle:(c,v)=>{if(v)this.classes.add(c);else this.classes.delete(c);}}; }
  appendChild(c){this.children.push(c);}
  addEventListener(k,fn){this.events[k]=fn;}
  showModal(){this.open=true;}
  close(){this.open=false;}
}
const elements = new Map();
function el(id){if(!elements.has(id))elements.set(id,new Element(id));return elements.get(id);}
const caseNodes=[];
const viewNodes=[new Element(),new Element()];viewNodes[0].dataset.view='top2';viewNodes[1].dataset.view='all';
const imageNodes=Array.from(html.matchAll(/data-image="([^"]+)"/g),m=>{const e=new Element();e.dataset.image=m[1];return e;});
const document={body:new Element('body'), getElementById:el,
  createElement:()=>{const e=new Element();caseNodes.push(e);return e;},
  querySelectorAll:q=>q==='.case'?caseNodes:q==='[data-view]'?viewNodes:q==='[data-image]'?imageNodes:[]};
const ctx = vm.createContext({document,console});
new vm.Script(script,{filename:'dashboard.html:script'}).runInContext(ctx);
const D=vm.runInContext('D',ctx), order=vm.runInContext('ORDER',ctx);
assert.strictEqual(document.body.dataset.currentCase,'competitive');
assert.strictEqual(document.body.dataset.representativeView,'top2');
assert.strictEqual(caseNodes.length,12);
assert.strictEqual(el('primaryButtons').children.length,6);
assert.strictEqual(el('refinedButtons').children.length,6);
assert.strictEqual((el('probRows').innerHTML.match(/class="prob-row"/g)||[]).length,2);
for(const p of ['50.6418%','40.4265%'])assert(el('probRows').innerHTML.includes(p));
assert(el('remainder').innerHTML.includes('8.9317%'));
const prim=JSON.parse(fs.readFileSync(path.join(root,'results/cases/competitive/probabilities.json'),'utf8'));
D.cases.competitive.p.forEach((p,i)=>assert(Math.abs(p-prim.result.topology_probabilities[i])<1e-15));
function pngInfo(b64){const b=Buffer.from(b64,'base64');assert(b.subarray(0,8).equals(Buffer.from([137,80,78,71,13,10,26,10])));return {width:b.readUInt32BE(16),height:b.readUInt32BE(20),bytes:b.length};}
let inlineImages=0;
for(const c of Object.values(D.cases)){
  for(const k of ['rep_top2','rep_all','unc']){pngInfo(c[k]);inlineImages++;}
  assert(Math.abs(c.p.reduce((a,b)=>a+b,0)-1)<1e-12);
  assert(Math.abs(c.remainder_modeled-c.p[c.order[2]])<1e-14);
}
for(const b64 of Object.values(D.images)){pngInfo(b64);inlineImages++;}
const widths=[];
for(const b of caseNodes){
  b.onclick();const c=D.cases[b.dataset.case];
  assert.strictEqual(document.body.dataset.currentCase,b.dataset.case);
  assert.strictEqual(el('repImg').src,'data:image/png;base64,'+c.rep_top2);
  assert.strictEqual(el('uncImg').src,'data:image/png;base64,'+c.unc);
  viewNodes[1].onclick();assert.strictEqual(document.body.dataset.representativeView,'all');
  assert.strictEqual(el('repImg').src,'data:image/png;base64,'+c.rep_all);
  assert(pngInfo(c.rep_top2).width<pngInfo(c.rep_all).width);
  widths.push({case:b.dataset.case,top2:pngInfo(c.rep_top2).width,all:pngInfo(c.rep_all).width});
  viewNodes[0].onclick();assert.strictEqual(document.body.dataset.representativeView,'top2');
}
caseNodes.find(b=>b.dataset.case==='competitive').onclick();
for(const b of imageNodes){b.onclick();assert.strictEqual(el('figureDialog').open,true);assert.strictEqual(el('zoomImg').src,el(b.dataset.image).src);el('closeZoom').onclick();assert.strictEqual(el('figureDialog').open,false);}
const result={status:'PASS_CODE_LEVEL',cases_checked:12,primary_cases:6,refined_cases:6,representative_modes_checked:2,
  inline_pngs_decoded:inlineImages,zoom_logic_checked:imageNodes.length,default_case:'competitive',default_top_k:2,renormalized:false,
  primary_percentages:{H1:'50.6418%',H2:'40.4265%',H0_remainder:'8.9317%'},remote_resources:0,
  browser_rendering:'NOT_RUN',browser_block:'CUA URL security policy rejected file:// and explicitly prohibited workaround paths',
  responsive_rules_present:true,responsive_visual_acceptance:'NOT_RUN',width_checks:widths,
  physics_calls:0,training_updates:0,
  wall_seconds:Number(process.hrtime.bigint()-started)/1e9,cpu_seconds:(process.cpuUsage(cpu0).user+process.cpuUsage(cpu0).system)/1e6};
fs.writeFileSync(path.join(__dirname,'html_logic_qa.json'),JSON.stringify(result,null,2)+'\n');
const rp=path.join(__dirname,'receipt.json');const rec=JSON.parse(fs.readFileSync(rp,'utf8'));rec.html_logic_qa={status:result.status,wall_seconds:result.wall_seconds,cpu_seconds:result.cpu_seconds,inline_pngs:result.inline_pngs_decoded,cases:12,browser_rendering:'NOT_RUN'};fs.writeFileSync(rp,JSON.stringify(rec,null,2)+'\n');
console.log(JSON.stringify({status:result.status,cases:12,images:inlineImages,defaultTop2:true,browser:'NOT_RUN',wall_seconds:result.wall_seconds,cpu_seconds:result.cpu_seconds}));
