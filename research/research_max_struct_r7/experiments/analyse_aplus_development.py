"""Describe paid disclosed A+ trajectories using public objectives only."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'research/research_max_struct_r7/delegated/aplus'

def main():
 results=json.loads((OUT/'DEVELOPMENT_RESULTS.json').read_text());rows=[]
 guarded=OUT/'GUARDED_DEVELOPMENT_RESULTS.json'
 if guarded.is_file():
  refinement=json.loads(guarded.read_text());results['rows']+=refinement['rows'];results['status']='COMPLETE' if results['status']==refinement['status']=='COMPLETE' else 'PARTIAL_REFINEMENT'
 arc=OUT/'ARC_DEVELOPMENT_RESULTS.json'
 if arc.is_file():
  proposal=json.loads(arc.read_text());results['rows']+=proposal['rows'];results['status']='COMPLETE' if results['status']==proposal['status']=='COMPLETE' else 'PARTIAL_ARC'
 for row in results['rows']:
  folder=Path(row['output_dir']);receipt=json.loads((folder/'RECEIPT.json').read_text());history=json.loads((folder/'HISTORY.json').read_text());iterations=[h for h in history if 'iteration'in h]
  records=[h for h in iterations if 'material_final_block_loss'in h]
  record=dict(case_id=row['case_id'],mode=row['mode'],full_cpu_seconds=row['cpu_seconds'],auxiliary_cpu_seconds=row['auxiliary_cpu_seconds'],original_Phi=row['original_Phi'],physical_relative_residual=row['physical_relative_residual'],physical_LU=row['physical_cost'].get('LU',0),physical_complex_RHS=row['physical_cost'].get('actual_complex_RHS',0),auxiliary_iterations=len(iterations),readout_loss_increases=sum(h['material_final_block_loss']>min(h['material_old_block_loss'],h['material_scalar_seed_block_loss'])+1e-14 for h in records),readout_fallbacks=sum(h['material_fallback'] for h in records),readout_calls=len(records),recorded_readout_cpu_seconds=sum(h.get('aplus_readout_cpu_seconds',0) for h in records),corrector=receipt.get('corrector'),selected_auxiliary_stage=receipt.get('auxiliary_selected_stage'),auxiliary_candidate_scores=receipt.get('auxiliary_candidate_scores'),output_dir=str(folder))
  rows.append(record)
 bykey={(r['case_id'],r['mode']):r for r in rows}
 for row in rows:
  base=bykey.get((row['case_id'],'original_a'))
  if base:
   row['Phi_over_original_A']=row['original_Phi']/base['original_Phi'];row['full_CPU_over_original_A']=row['full_cpu_seconds']/base['full_cpu_seconds'];row['auxiliary_CPU_over_original_A']=row['auxiliary_cpu_seconds']/base['auxiliary_cpu_seconds']
 (OUT/'ANALYSIS.json').write_text(json.dumps(dict(status=results['status'],truth_access='NONE',scientific_scope='Original physical Phi and cost retention only; GT recovery NOT_SCORED; disclosed development is not blind confirmation.',rows=rows),indent=2)+'\n')
 lines=['# A+ disclosed full-trajectory development','',f"Status: {results['status']}; truth access NONE; recovery scoring NOT_RUN.",'','All runs use the same 3000 auxiliary iterations and original 3000-iteration corrector allowance. Ratios compare to the same-object original-A control. Original Phi is the online observation/prior objective, not a GT recovery metric.','','| Object | Readout | CPU s | Aux CPU s | Phi | Phi / A | CPU / A | LU | RHS | inner loss increases |','|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
 for r in rows:
  lines.append(f"| {r['case_id']} | {r['mode']} | {r['full_cpu_seconds']:.3f} | {r['auxiliary_cpu_seconds']:.3f} | {r['original_Phi']:.8g} | {r.get('Phi_over_original_A',1):.6f} | {r.get('full_CPU_over_original_A',1):.3f} | {r['physical_LU']} | {r['physical_complex_RHS']} | {r['readout_loss_increases']} |")
 lines+=['','Each output directory contains full auxiliary scalar histories, stage joint m/J/Q/A, paid physical stage states and selected/final predictions, corrector evaluation/accepted histories, and process CPU/LU/RHS receipts. Source and configuration are frozen separately.','']
 (OUT/'FULL_TRAJECTORY_COMPARISON.md').write_text('\n'.join(lines))
 print('\n'.join(lines))
if __name__=='__main__':main()
