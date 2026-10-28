"""Build small reading arrays only AFTER prediction seal and offline scoring."""
from pathlib import Path
import json,hashlib,numpy as np
P=Path(__file__).resolve().parents[1]
def main():
    assert (P/'E2/PREDICTIONS_SEALED.json').is_file() and (P/'E2/OFFLINE_GT_SCORES.json').is_file()
    parents=[f'R7_LOCKED_{i:03d}' for i in range(8)]
    arms=['A_FROZEN','B_FULLJ','C_JETS','D_RANDOM']
    truth=[];obs=[];anchors=[];xx=[];JJ=[];mm=np.full((4,8,576),np.nan+1j*np.nan,complex);pred=np.full((4,8,32,8),np.nan+1j*np.nan,complex)
    available=np.zeros((4,8),bool);phi=np.full((4,8),np.nan);checks=0
    for i,parent in enumerate(parents):
        d=P/'E2/parents'/parent
        with np.load(P/'raw/offline_truth'/f'{parent}.npz',allow_pickle=False) as z:truth.append(z['material'].copy())
        with np.load(d/'OBSERVATION.npz',allow_pickle=False) as z:obs.append(z['observation'].copy())
        with np.load(d/'ANCHOR.npz',allow_pickle=False) as z:anchors.append(z['material'].copy());xx.append(z['x'].copy());JJ.append(z['J'].copy())
        for j,arm in enumerate(arms):
            f=d/arm/'RMS_0.0200/ENDPOINT.npz'
            if not f.exists():continue
            with np.load(f,allow_pickle=False) as z:mm[j,i]=z['material'];pred[j,i]=z['prediction'];available[j,i]=True
            r=json.loads((f.parent/'RESULT.json').read_text());phi[j,i]=r['Phi']
    a=dict(parents=np.array(parents),arms=np.array(arms),GT_material=np.array(truth),observation=np.array(obs),anchor_material=np.array(anchors),anchor_x=np.array(xx),anchor_J=np.array(JJ),candidate_material=mm,candidate_prediction=pred,available=available,Phi=phi)
    out=P/'GPT_PRO_HANDOFF/compact';out.mkdir(exist_ok=True)
    path=out/'R7_LOCKED_COMPACT.npz';np.savez_compressed(path,**a)
    with np.load(path,allow_pickle=False) as z:
        for k,v in a.items():
            assert z[k].dtype==v.dtype and z[k].shape==v.shape
            if v.dtype.kind in 'fc':assert np.array_equal(z[k],v,equal_nan=True)
            else:assert np.array_equal(z[k],v)
            checks+=1
    (out/'MANIFEST.json').write_text(json.dumps(dict(scope='exposed post-seal compact reading arrays, NOT new blind data or full raw evidence',bytes=path.stat().st_size,blake2b=hashlib.blake2b(path.read_bytes()).hexdigest(),array_readbacks=checks,available_counts=available.sum(1).tolist(),missing='C003/004 use NaN with available=false, never copied anchor padding',GT_inclusion='offline reading and scoring only'),indent=2)+'\n')
    print('COMPACT_EXACT_READBACK_PASS',checks,path.stat().st_size)
if __name__=='__main__':main()
