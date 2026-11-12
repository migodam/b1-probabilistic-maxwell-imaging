"""Retrospective oracle diagnostic: did AP ever *generate* a good audited material?
Truth is used only after frozen predictions for analysis, never to change outputs.
"""
from pathlib import Path
import os
import json
import numpy as np
import pandas as pd
ROOT=Path(os.environ.get('MAX_STRUCT_RUN_ROOT', str(Path(__file__).resolve().parent.parent)))
BASE=Path(os.environ.get('MAX_STRUCT_REPO', str(ROOT/'research_repo')))/'research/research_max_struct_r12/supplement_paired'
INPUT=ROOT/'work/paired_selected_candidates'
OUT=ROOT/'work/evidence_analysis'
rows=[];traces=[]
for path in sorted(INPUT.glob('*__CANDIDATES.npz')):
    task=path.name.removesuffix('__CANDIDATES.npz')
    cohort,case,method=task.split('__')
    history=pd.DataFrame(json.loads((BASE/'RESULTS'/task/'HISTORY.json').read_text()))
    with np.load(path) as z:M=z['m']
    assert len(M)==len(history)
    with np.load(BASE/'offline_truth'/cohort/f'{case}.npz') as z: gt=z['m'];amp=float(z['amplitude'])
    history['offline_NRMSE']=np.linalg.norm(M-gt[None,:],axis=1)/(np.sqrt(len(gt))*amp)
    history['cohort']=cohort;history['case_id']=case;history['method']=method
    eligible=history[history.eligible]
    proposal=eligible[eligible.phase.isin(['initial','proposal'])]
    selected=proposal.loc[proposal.Phi.idxmin()]
    oracle=proposal.loc[proposal.offline_NRMSE.idxmin()]
    alloracle=eligible.loc[eligible.offline_NRMSE.idxmin()]
    final=eligible.loc[eligible.Phi.idxmin()]
    row=dict(task_id=task,cohort=cohort,case_id=case,method=method,
      proposal_candidates=len(proposal),selected_raw_nrmse=selected.offline_NRMSE,
      oracle_proposal_nrmse=oracle.offline_NRMSE,oracle_proposal_candidate=int(oracle.candidate),
      oracle_proposal_CPU_seconds=oracle.CPU_seconds,oracle_proposal_Phi=oracle.Phi,
      oracle_proposal_residual_over_delta=oracle.residual_over_delta,
      oracle_all_nrmse=alloracle.offline_NRMSE,primary_nrmse=final.offline_NRMSE,
      selected_raw_candidate=int(selected.candidate))
    rows.append(row)
    traces.append(history[['cohort','case_id','method','candidate','CPU_seconds','phase','Phi','D','R','residual_over_delta','eligible','offline_NRMSE']])
pd.DataFrame(rows).to_csv(OUT/'candidate_access_oracle_diagnostic.csv',index=False)
pd.concat(traces,ignore_index=True).to_csv(OUT/'candidate_access_traces_offline.csv',index=False)
print(pd.DataFrame(rows)[['case_id','method','proposal_candidates','selected_raw_nrmse','oracle_proposal_nrmse','oracle_proposal_residual_over_delta','oracle_all_nrmse','primary_nrmse']].to_string(index=False))
