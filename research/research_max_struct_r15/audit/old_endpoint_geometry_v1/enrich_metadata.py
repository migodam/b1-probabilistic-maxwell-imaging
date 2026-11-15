"""Metadata/accounting completion only; no array/physics recomputation."""
from pathlib import Path
import json,csv,zipfile,struct,ast,time
from collections import Counter
out=Path(__file__).resolve().parent
cpu,wall=time.process_time(),time.perf_counter()
def write(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def tab(p,rows):
 fields=list(dict.fromkeys(k for r in rows for k in r))
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
  for r in rows:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in r.items()})
shapes={}
for p in sorted((out/'public_inputs').glob('*.npz')):
 vals={}
 with zipfile.ZipFile(p) as z:
  for n in ['E0','S']:
   with z.open(n+'.npy') as f:
    f.read(6);version=f.read(2);size=struct.unpack('<H' if version[0]==1 else '<I',f.read(2 if version[0]==1 else 4))[0]
    vals[n]=ast.literal_eval(f.read(size).decode('latin1'))['shape']
 shapes[p.stem]=dict(sources=vals['E0'][1],receivers=vals['S'][0],material_dof=vals['E0'][0])
rows=json.loads((out/'ENDPOINTS.json').read_text());total=Counter();percase={}
for r in rows:
 dimensions=shapes[r['case_id']];r.update(dimensions);ns=dimensions['sources']
 extra=dict(OFFLINE_dense_G_output_matrix_products=2,OFFLINE_dense_G_output_complex_columns=8*ns,
            OFFLINE_dense_S_output_matrix_products=4,OFFLINE_dense_S_output_complex_columns=10*ns,
            OFFLINE_dense_S_material_null_matrix_products=1,OFFLINE_dense_S_material_null_columns=1,
            OFFLINE_saved_basis_orthogonal_projection_actions=14)
 r['physical_cost'].update(extra)
 percase.setdefault(r['case_id'],Counter()).update(extra);total.update(extra)
 write(out/'cases'/r['case_id']/(r['method']+'_RECEIPT.json'),r)
write(out/'ENDPOINTS.json',rows);tab(out/'ENDPOINTS.csv',rows)
casecost=json.loads((out/'CASE_COSTS.json').read_text())
for r in casecost:r['physical_cost'].update(percase[r['case_id']])
write(out/'CASE_COSTS.json',casecost)
for case,d in shapes.items():
 p=out/'cases'/case/'CASE_METADATA.json';v=json.loads(p.read_text());v.update(d);write(p,v)
# Regenerate event table from actual preserved rows; shared GT costs stay unchanged.
with (out/'COST_EVENTS.csv').open() as f:events=list(csv.DictReader(f))
by={(r['case_id'],r['method']):r for r in rows}
for e in events:
 if (e['case_id'],e['method']) in by:e['cost']=by[e['case_id'],e['method']]['physical_cost']
 elif isinstance(e['cost'],str):e['cost']=json.loads(e['cost'])
tab(out/'COST_EVENTS.csv',events)
p=out/'RECEIPT.json';receipt=json.loads(p.read_text());receipt['physical_cost'].update(total)
receipt['public_sources_per_case']=shapes
receipt['source_count_distribution']={str(k):v for k,v in Counter(d['sources'] for d in shapes.values()).items()}
receipt['supplementary_counter_scope']='Dense output products and saved-basis projection counts are exact source-level counts of the already executed computation, enriched without rerunning arrays/physics. Matrix-product column counts are not LU-solved actual_complex_RHS.'
write(p,receipt)
p=out/'REPORT.md';text=p.read_text();text+='\nActual source counts are heterogeneous:16 cases have12 sources,8 have6. All have48 receivers. Costs use actual ns.72LU and6000 LU-solved complexRHS include shared GT once/case and separate exact-Snull actions. Additional dense G/S output products and saved-basis projection counts are listed separately from LU RHS in RECEIPT.\n'
p.write_text(text)
write(out/'METADATA_ENRICHMENT_RECEIPT.json',dict(status='COMPLETE',science_rerun=False,array_payloads_read=False,source='Saved public NPZ key headers and exact executed-source product counts',CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall))
print('metadata enriched; no numerical recomputation')
