"""Post-seal explanatory identity; no model or primary hypothesis changes."""
from pathlib import Path
import json,time
import numpy as np
P=Path(__file__).resolve().parents[1]
def main():
    if not (P/'raw/predictions/ALL_PREDICTIONS_SEALED.json').exists():raise RuntimeError('sealed images required')
    with np.load(P/'raw/offline_scores/OFFLINE_SCOREPACK.npz',allow_pickle=False) as z:
        methods=z['methods'].astype(str).tolist();truth=z['truth_chi'];raw=z['prediction_raw'];clip=z['prediction_clipped']
    result=dict(role='EXPLORATORY_AFTER_ALL_PREDICTIONS_FROZEN',identity='mean_s||m_s-truth||^2=||mean_s m_s-truth||^2+mean_s||m_s-mean_s m_s||^2',
                loss='complex full-grid /576',seeds=3,parents=512,utc=time.time(),no_new_training=True,new_RHS=0,
                caveat='finite three-checkpoint dispersion; not a population training-variance estimate or an alternative primary retention test')
    for label,preds in [('raw',raw),('clipped',clip)]:
        mean=preds.mean(1);avg=np.abs(preds-truth[:,None,None,:])**2;ens=np.abs(mean-truth[:,None,:])**2
        dispersion=(np.abs(preds-mean[:,None,:,:])**2).mean((1,3))
        avg=avg.mean((1,3));ens=ens.mean(2)
        error=np.max(np.abs(avg-ens-dispersion));assert error<1e-15
        means={m:dict(seed_mean_risk=float(avg[:,j].mean()),ensemble_risk=float(ens[:,j].mean()),three_checkpoint_dispersion=float(dispersion[:,j].mean())) for j,m in enumerate(methods)}
        k,d=methods.index('k32'),methods.index('distilled_w128')
        gap=float((avg[:,d]-avg[:,k]).mean());eg=float((ens[:,d]-ens[:,k]).mean());vg=float((dispersion[:,d]-dispersion[:,k]).mean())
        result[label]=dict(methods=means,max_identity_error=float(error),distilled_minus_k32=dict(seed_mean_gap=gap,ensemble_gap=eg,dispersion_gap=vg,fraction_gap_due_to_dispersion=vg/gap if gap else None))
    (P/'ENSEMBLE_VARIANCE_DECOMPOSITION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['raw']['distilled_minus_k32']))
if __name__=='__main__':main()
