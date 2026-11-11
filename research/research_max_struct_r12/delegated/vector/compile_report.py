"""Post-seal vector figures/readout; no forward solve or optimization."""
from pathlib import Path
import json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).absolute().parent
data=json.loads((ROOT/'SCORES.json').read_text());rows=data['rows'];cases=list(dict.fromkeys(r['case_id'] for r in rows));methods=list(dict.fromkeys(r['method'] for r in rows))
longrows=json.loads((ROOT/'LONG_SCORES.json').read_text())['rows']
fig,axs=plt.subplots(len(cases),len(methods)+1,figsize=(18,8),constrained_layout=True)
for i,cid in enumerate(cases):
    truth=np.load(ROOT/'sealed_truth'/(cid+'.npz'))['chi']
    for j,method in enumerate(['truth']+methods):
        m=truth if method=='truth' else np.load(ROOT/'RESULTS'/cid/method/'prediction.npz')['prediction'].real
        ax=axs[i,j];im=ax.imshow(m.reshape(6,6,6)[:,:,3].T,origin='lower',vmin=0,vmax=4,cmap='magma',extent=[-.3,.3,-.3,.3]);ax.set_xticks([]);ax.set_yticks([])
        if i==0:ax.set_title(method)
        if j==0:ax.set_ylabel(cid.replace('_','\n'),fontsize=8)
        if method!='truth':
            r=next(r for r in rows if r['case_id']==cid and r['method']==method);ax.set_xlabel(f"NRMSE {r['NRMSE']:.3f}\nr={r['physical_data_residual']:.4f}",fontsize=8)
fig.colorbar(im,ax=axs,label='real material contrast chi; mid-z slice',shrink=.6)
fig.suptitle('Reduced true3D: 216 unknown material cells, 648 vector field components; offline scores',fontsize=12)
fig.savefig(ROOT/'MATERIAL_SLICES.png',dpi=170);plt.close(fig)
fig,axs=plt.subplots(len(cases),len(methods)+1,figsize=(18,8),constrained_layout=True)
for i,cid in enumerate(cases):
    truth=np.load(ROOT/'sealed_truth'/(cid+'.npz'))['chi']
    for j,method in enumerate(['truth']+methods):
        m=truth if method=='truth' else np.load(ROOT/'RESULTS'/cid/method/'long_corrector/prediction.npz')['prediction']
        ax=axs[i,j];im=ax.imshow(m.reshape(6,6,6)[:,:,3].T,origin='lower',vmin=0,vmax=4,cmap='magma',extent=[-.3,.3,-.3,.3]);ax.set_xticks([]);ax.set_yticks([])
        if i==0:ax.set_title(method)
        if j==0:ax.set_ylabel(cid.replace('_','\n'),fontsize=8)
        if method!='truth':
            r=next(r for r in longrows if r['case_id']==cid and r['method']==method);ax.set_xlabel(f"NRMSE {r['NRMSE']:.3f}\nr={r['physical_data_residual']:.4f}",fontsize=8)
fig.colorbar(im,ax=axs,label='real material contrast chi; mid-z slice',shrink=.6)
fig.suptitle('Same six methods after identical longer physical continuation (max1000 iterations /60CPU s)',fontsize=12)
fig.savefig(ROOT/'MATERIAL_SLICES_LONG.png',dpi=170);plt.close(fig)
costrows=[]
for r in rows:
    rec=json.loads((ROOT/'RESULTS'/r['case_id']/r['method']/'RECEIPT.json').read_text());row=dict(case_id=r['case_id'],method=r['method'],CPU_seconds=rec['whole_child_CPU_seconds_near_exit'],physical_counts=rec['counts'])
    if r['method'] in ['AP','DR']:
        hist=json.loads((ROOT/'RESULTS'/r['case_id']/r['method']/'projection/history.json').read_text());row.update(projection_setup_cpu=rec['projection_setup_cpu'],projection_cpu=rec['projection_cpu'],projection_counts=rec['projection_counts'],projection_setup_counts={k:hist['history'][0]['costs'][k] for k in ['H_cholesky','H_setup_rhs','port_eigh']},iterations=rec['projection_iterations'],archive_and_corrector_CPU_residual=rec['whole_child_CPU_seconds_near_exit']-rec['projection_setup_cpu']-rec['projection_cpu'])
    costrows.append(row)
generation=json.loads((ROOT/'GENERATION_RECEIPT.json').read_text())
dump=dict(status='COMPLETE',solver_CPU_seconds_lower_bound=sum(r['CPU_seconds'] for r in costrows),generation_CPU_seconds=generation['CPU_seconds'],rows=costrows,scope='Per-arm timing starts before model setup and includes state archives and final outputs except final receipt write; initial module imports excluded in first runs. This sum is a lower bound. Queue external child CPU where present is separate not additive.',work_balance='3000 AP/DR projections vs120 inherited PR cycles; fixed declared work, not equal CPU; identical120-maxiter common corrector. Physical LU counts exclude fixed-background H Cholesky and source/receiver products.')
dump['long_continuation_CPU_seconds_lower_bound']=sum(r['additional_CPU_seconds'] for r in longrows)
dump['long_continuation_counts']=[dict(case_id=r['case_id'],method=r['method'],counts=r['counts']) for r in longrows]
dump['total_CPU_seconds_lower_bound']=dump['solver_CPU_seconds_lower_bound']+dump['long_continuation_CPU_seconds_lower_bound']+generation['CPU_seconds']
(ROOT/'COST_LEDGER.json').write_text(json.dumps(dump,indent=2)+'\n')
lines=['# Reduced true3D exploratory material inversion','', 'Three public development observations: matched mixed coated cube, matched connected unequal inclusions, and the same coated cube generated on an8³ grid and inverted on6³. All inverse grids have216 real chi unknowns and648 complex vector-current components, with60 common dipole sources and180 electric receiver channels. Bounds chi[0,4] (epsilon[1,5]), uniform chi=.25 and1% raw iid circular-complex noise are fixed. This is a reduced dyadic CM+RR DDA adaptation, not a literal2019 reproduction.','', 'Baseline methods directly call inherited R8 VectorPR3D.optimize_pr: CIE_sourcewise, Dual_sourcewise, CIE_global and Aplus_global. Their120 auxiliary cycles use inherited mf1→2 and beta6→1 (Dual current beta3 late), inner12. AP/DR use scale6, 3000 projections/audit50/320CPU cap with one raw iid99% noise ball and no whitening. Physical candidate selection and all common correctors use the same inherited source-normalized physical Phi and same quadratic prior, with corrector maxiter120. The projection feasibility ball and candidate-scoring objective have explicitly different data normalization.','', 'All18 predictions were sealed before scorer GT read. Selection never uses GT. Crossgrid truth for scoring is the same analytic geometry averaged onto the inverse cells. No independent fresh confirmation or high-epsilon8 claim.','', '| Case | Method | Final NRMSE | Raw NRMSE | Physical residual | Phi | Paid CPU s |','|---|---|---:|---:|---:|---:|---:|']
for r in rows:lines.append(f"| {r['case_id']} | {r['method']} | {r['NRMSE']:.5f} | {r['raw_NRMSE']:.5f} | {r['physical_data_residual']:.5f} | {r['physical_phi']:.6g} | {r['CPU_seconds']:.2f} |")
lines += ['', 'All original120 common correctors hit their iteration cap. Before any GT score, root requested an identical additional max1000iterations/60CPU seconds physical continuation for all18 endpoints; CONTINUATION_CONTRACT.json records that common budget diagnostic. Original120 outputs and costs remain intact. Both18+18 endpoints were sealed before scoring.','', '| Case | Method | Longer-corrector NRMSE | Physical residual | Phi | Extra CPU s | Stop |','|---|---|---:|---:|---:|---:|---|']
for r in longrows:lines.append(f"| {r['case_id']} | {r['method']} | {r['NRMSE']:.5f} | {r['physical_data_residual']:.5f} | {r['physical_phi']:.7g} | {r['additional_CPU_seconds']:.2f} | {r['optimizer']['message']} |")
lines += ['', 'A relative-function optimizer stop is not a local-minimum certificate. The original120 stage and longer common continuation answer different budget questions. Neither extending the common corrector nor a lower objective establishes superior material recovery by itself.']
lines += ['', 'CPU is measured work, not equal-compute evidence. All fixed-background Cholesky/setup/eigensystems, actions, projection multiplier evaluations, physical audits, common correctors and saved full audit-state archives are recorded separately in COST_LEDGER.json and receipts. The declared work differs substantially across method families; paid prefix histories permit physical-Phi comparisons but no claim of equal CPU superiority is made. Archive serialization is charged and can be substantial. Numerical recovery/residual and scientific mechanism support remain separate judgments.','', 'Evidence: CONTRACT.json, PROJECTION_CONTRACT.json, GENERATION_RECEIPT.json, data_public/, sealed_truth/, RESULTS/{case}/{method}/prediction.npz and RECEIPT.json, projection/history.json and states.npz, PREDICTION_SEAL.json, SCORES.json, COST_LEDGER.json, MATERIAL_SLICES.png. No new physical solve occurs in this report compiler.']
(ROOT/'REPORT.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(dict(status='COMPLETE',rows=len(rows),solver_CPU_lower_bound=dump['solver_CPU_seconds_lower_bound'])))
