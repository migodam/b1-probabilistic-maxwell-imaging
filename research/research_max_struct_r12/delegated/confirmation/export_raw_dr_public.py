"""Extract saved DR audit materials only; no GT, projection, physics or selection tuning."""
from pathlib import Path
import json,zipfile,io,time
import numpy as np
ROOT=Path(__file__).resolve().parent
def main():
    target=ROOT/'R12_RAW_DR_AUDIT_MATERIALS_PUBLIC.zip'
    if target.exists():raise FileExistsError(target)
    index=[]
    with zipfile.ZipFile(target,'x',compression=zipfile.ZIP_DEFLATED) as archive:
        for stage in ('dev6','confirm24'):
            seal=json.loads((ROOT/('SEALED_'+stage.upper()+'_PUBLIC_COHORT.json')).read_text())
            for row in seal['rows']:
                if row['method']!='dr':continue
                task=ROOT/'RESULTS'/stage/row['task_id']
                core=json.loads((task/'projection/history.json').read_text())
                physical=json.loads((task/'history.json').read_text())
                # Search runner performs one initial audit before constructing
                # the projection, then projection audit0 and scheduled audits.
                recs=[];materials=[]
                with np.load(task/'projection/states.npz',allow_pickle=False) as z:
                    for k,c in enumerate(core['history']):
                        p=physical[k+1]
                        if float(c['objective'])!=float(p['Phi']):raise ValueError('Audit ordering mismatch')
                        materials.append(z[f'audit_{k:05d}_candidate_material'].copy())
                        recs.append(dict(iteration=c['iteration'],global_CPU_seconds=p['CPU_seconds'],
                            Phi=p['Phi'],R=p['prior'],relative_data_residual=p['relative_data_residual'],
                            noise_feasible=p['noise_feasible'],eligible_proposal144=bool(p['eligible'] and p['CPU_seconds']<=144.),
                            gap=c['gap'],source_audit_index=k))
                dest=Path(stage)/row['task_id'];buf=io.BytesIO()
                np.savez_compressed(buf,candidate_materials=np.asarray(materials))
                archive.writestr(str(dest/'raw_dr_audits.npz'),buf.getvalue())
                data=dict(case_id=row['case_id'],stage=stage,method='dr',records=recs,
                    iterations_completed=core['iterations'],stop_reason=core['stop_reason'],
                    source=str(task/'projection/states.npz'),new_FW=0,GT='NONE',
                    scope='Saved physical-audited C-shadow material; no next-step/un-audited shadow reconstructed')
                archive.writestr(str(dest/'raw_dr_audits.json'),json.dumps(data,indent=2))
                index.append(dict(case_id=row['case_id'],stage=stage,npz=str(dest/'raw_dr_audits.npz'),metadata=str(dest/'raw_dr_audits.json')))
        archive.writestr('RAW_DR_PUBLIC_INDEX.json',json.dumps(dict(rows=index,new_FW=0,GT='NONE',process_CPU_seconds_near_write=time.process_time()),indent=2))
    print(json.dumps(dict(status='EXTRACTED_EXISTING_AUDITS_ONLY',cases=len(index),bytes=target.stat().st_size,CPU=time.process_time(),new_FW=0,GT='NONE')),flush=True)
if __name__=='__main__':main()
