"""Small independent saved-array replay. No PDE, truth-based model selection or hashes."""
from pathlib import Path
import json,re,ast,csv,numpy as np
from r2.accounting import BudgetBook
ROOT=Path(__file__).resolve().parents[1]
def main():
    required=['START_HERE.md','CONTRACT.json','PRE_REGISTERED_PLAN.md','THEORY_TO_TEST.md','EVIDENCE_LEDGER.csv','QUERY_LEDGER.jsonl','PHYSICS_COST.json','E1_RESULTS.csv','E2_RESULTS.csv','E3_RESULTS.csv','THEOREM_STATUS.md','GO_NO_GO.json','REPORT_ZH.md','NEXT_PRO_HANDOFF.md','REPRODUCE.md']
    assert all((ROOT/p).is_file() for p in required)
    parse=list((ROOT/'src').rglob('*.py'))+list((ROOT/'tools').rglob('*.py'))
    for p in parse:ast.parse(p.read_text())
    micro=[]
    for parent in [0,1]:
        a=np.load(ROOT/f'results/micro/WALL_p{parent}_n16_t0.npz');b=np.load(ROOT/f'results/micro/WALL_p{parent}_n16_t1.npz');mask=a['eps']!=1
        rms=float(np.sqrt(np.mean((a['eps'][mask]-b['eps'][mask])**2)));assert abs(rms-np.sqrt(1.28))<1e-12
        for j in range(12):
            arr=[]
            for tex in [0,1]:
                fs=[]
                for f in [.85,1.,1.15]:
                    p=ROOT/f'results/micro/p{parent}_n16_t{tex}_f{f:g}_q{j+2}.npz'
                    if not p.exists():p=ROOT/f'results/imaging/micro_p{parent}_t{tex}_f{f:g}_q{j}.npz'
                    fs.append(np.load(p)['data'])
                arr.append(np.asarray(fs))
            micro.append(float(np.linalg.norm(arr[0]-arr[1])/.01))
    f=ROOT/'results/final_pair_refinement';d=[]
    for hz in [.85,1.,1.15]:
        p=np.load(f/f'f{hz}.npz');d.append(p['endpoint']-p['base'])
    pair=float(np.linalg.norm(d)/.01);saved=json.loads((f/'RESULTS.json').read_text())['joint_distance'];assert abs(pair-saved)<1e-12
    b=BudgetBook();caps=b.auth
    for phase,c in caps['caps'].items():
        x=b.snapshot(phase);assert x['rhs']<=c['rhs'] and x['full_fe']<=c['full_fe'] and x['process_cpu_seconds']<=c['process_cpu_seconds']
    links=[]
    for p in ROOT.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if link.startswith(('http:','https:','#')):continue
            q=link.split('#')[0]
            if q and not (p.parent/q).exists():links.append(dict(file=p.name,link=link))
    assert not links,links
    result=dict(status='PASS_SAVED_REPLAY_NOT_NEW_MAXWELL_QUALIFICATION',required_files=len(required),parsed_source_files=len(parse),heldout_fine_A2_max_white=max(micro),material_RMS=np.sqrt(1.28),final_pair_recomputed=pair,final_pair_saved=saved,source_parse=True,missing_links=links,new_FE=0,new_RHS=0,no_hashes=True,PNG_count=len(list((ROOT/'figures').glob('*.png'))),PDF_count=len(list((ROOT/'figures').glob('*.pdf'))))
    (ROOT/'results/SAVED_REPLAY.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
if __name__=='__main__':main()
