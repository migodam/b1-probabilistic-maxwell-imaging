"""Three matched 64-call corrected proposals, plus scoped cold/warm timing."""
import os
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
C1=Path(os.environ.get('C1_P2',ROOT.parent.parent/'C1/c1-joint-p2'));sys.path.insert(0,str(C1/'p2'))
from src.b1_joint.accounting import Ledger,dump
from src.b1_joint.forward import MaterialForward
from src.b1_joint.proposal_correction import run_proposal_correction,time_proposals
def main():
    out=ROOT/'research/joint_p2/results';ledger=Ledger(out);ledger.stage='reserve';f=MaterialForward(ROOT,ledger);models=C1/'p2/results/models';item=next(a for a in json.loads((out/'DEVELOPMENT_CASES.json').read_text())['cases'] if a['id']=='large_cavity');y=np.array(item['y'])[:32];sigma=np.array(item['sigma'])[:32]
    with ledger.job('p2-three-matched-full-likelihood-proposals'):
        summary,arrays=run_proposal_correction(f,models,y,sigma,json.loads((out/'FULL_MIE_CORRECTED_REFERENCE.json').read_text()))
        for method,a in arrays.items():np.savez_compressed(out/f'CORRECTED_PROPOSAL_{method}.npz',**a)
        dump(out/'CORRECTED_PROPOSALS.json',summary)
        dump(out/'PROPOSAL_TIMING.json',time_proposals(models,y,sigma))
    dump(out/'RESOURCE_SUMMARY.json',ledger.summary())
if __name__=='__main__':main()
