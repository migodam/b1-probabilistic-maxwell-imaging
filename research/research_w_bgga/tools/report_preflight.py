"""Publishable factual readiness report and geometry figure; no scientific solve."""
from pathlib import Path
import ast,csv,json,re
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bgga.accounting import BudgetBook

ROOT=Path(__file__).resolve().parents[1]
def dump(name,x):(ROOT/name).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')

def main():
    if any(r.get('rhs',0) or r.get('factorizations',0) for r in BudgetBook().records()):
        raise RuntimeError('Do not overwrite scientific phase outputs with a preflight-only report')
    b=json.loads((ROOT/'results/preflight/MESH_BENCHMARK.json').read_text())
    with np.load(ROOT/'results/preflight/FIXED_MESH.npz',allow_pickle=False) as z:
        p,t,q=z['points'],z['triangles'],z['target_quad_points']
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    axes[0].triplot(p[0],p[1],t.T,lw=.15,color='#a2afb8')
    axes[0].set(xlim=(-2,2),ylim=(.2,3),xlabel='x / wavelength',ylabel='y / wavelength',title='Actual public mesh; no target boundary fitted')
    axes[0].add_patch(plt.Rectangle((-1.5,.5),3,.75,fill=False,edgecolor='#c45842',lw=1.4,label='Public wall ROI'))
    axes[0].add_patch(plt.Rectangle((-1.2,1.4),2.4,1.4,fill=False,edgecolor='#166e87',lw=1.5,label='Entire target ROI'))
    axes[0].scatter(q[::20,0],q[::20,1],s=1,color='#1b9668',label='Example cut-polygon quadrature')
    axes[0].set_aspect('equal');axes[0].legend(fontsize=7)
    axes[1].bar(['Whole target ROI','Example disk'],[b['public_ports'],b['example_touched_ports']],color=['#166e87','#1b9668'])
    axes[1].set(ylabel='Structural FE coordinates (count)',title='Port dimension is not material parameter count')
    for i,v in enumerate([b['public_ports'],b['example_touched_ports']]):axes[1].text(i,v+25,str(v),ha='center')
    axes[1].text(.5,.52,'Full Green actions: 1.03 GB\nPublic setup: 2257 RHS\nDirect new target: 2 source RHS\n\nCounts only; accuracy and\ncold-time savings NOT TESTED',ha='center',transform=axes[1].transAxes,fontsize=10)
    fig.savefig(ROOT/'figures/PUBLIC_MESH_AND_PORTS.png',dpi=170);fig.savefig(ROOT/'figures/PUBLIC_MESH_AND_PORTS.pdf');plt.close(fig)
    dump('figures/PUBLIC_MESH_AND_PORTS_DATA.json',b)
    events=BudgetBook().records()
    measured=sum(x.get('process_cpu_seconds',0.) for x in events)
    dump('COST_SUMMARY.json',dict(status='PREFLIGHT_ONLY',new_full_fe=sum(x.get('factorizations',0) for x in events),
        new_complex_rhs=sum(x.get('rhs',0) for x in events),measured_wrapper_process_cpu_seconds=measured,
        authority='PHYSICS_LEDGER.jsonl job_completion; each once, kernel timers not added',
        wrapper_wall_seconds=sum(x.get('wall_seconds',0.) for x in events if x.get('label')=='job_completion'),
        peak_measured_job_RSS_bytes=max([x.get('children_peak_RSS',0) for x in events]+[0]),
        pending_current_report_job='its final receipt supersedes this in-child snapshot',
        accounting_gap='Early source/PDF extraction and worker static reads were not wrapped; not zero CPU. No physical call occurred.',
        standalone_deployment_cost='NOT_MEASURED',physical_speed_claim=False,
        historical_bank=dict(rhs=11800,fe=2950,cpu_seconds=9288.02304,charged_to_new_campaign=False)))
    dump('GATE_DECISION.json',dict(overall='AWAITING_G0_RESOURCE_AUTHORIZATION',
        EXACT_DISCRETE_ASSIMILATION='NOT_RUN_MAXWELL; synthetic matrix identity tested',
        TARGET_INVARIANT_PROPAGATION_CACHE='IMPLEMENTED_NOT_PHYSICALLY_QUALIFIED',
        ANGLE_ONLY_SAFETY='not generally sufficient analytically; new 12-parent Maxwell panel NOT_RUN',
        TASK_CONDITIONAL_PROPAGATION_IDENTIFIABILITY='HISTORICAL_CONDITIONAL_DISCRETE_EVIDENCE_ONLY',
        FINITE_PRIOR_POSTERIOR='NOT_RUN; interval algebra verified on small finite tests only',
        JOINT_DEPENDENCE='implemented preservation/lift; 2950-atom method comparison NOT_RUN',
        A0_UNKNOWN_SUPPORT='NOT_RUN_NEW; inherited A0_ONLINE_NOT_REACHED unchanged',
        COLD_COST='NOT_ESTABLISHED',CONTINUUM_CERTIFICATE='NOT_ESTABLISHED',
        resource_gate='No new BGGA grant in sources; user prompt explicitly requires it',
        scientific_negative_not_inferred_from_resource_stop=True))
    dump('FACTOR_HISTORY.json',dict(scope='STATIC_BGGA_V1',actual_observation_factors=[],
        status='NOT_RUN',policy='same data factor is version-replaced; old and new are never multiplied together',
        G1_planned=['obs_base_tx0_all16_rx','obs_base_tx1_all16_rx'],
        current_data_static=True,acquisition_panel='NOT_RUN'))
    with (ROOT/'FULL_TABLE.csv').open('w',newline='') as f:
        fields=['stage','parent','noise','method','seed','status','reason']
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        for stage in ['G0','G1']:
            methods=['FIXED_FULL','INTERFACE_FITTED','EXACT_RESOLVENT','EXACT_SCHUR','ORDINARY_KRYLOV','ORDINARY_GOAL_ROM'] if stage=='G0' else ['M0_DIRECT','M1_CONFIDENCE','M2_RETENTION','M3_TASK_POINT','M4_MULTI_GREEN','M5_QUOTIENT','M6_ORDINARY_CACHE_ROM']
            for parent in range(12):
                for noise in ([0] if stage=='G0' else [0,1]):
                    for method in methods:writer.writerow(dict(stage=stage,parent=parent,noise=noise,method=method,seed='predeclared',status='NOT_RUN_RESOURCE',reason='new physics grant absent; no fabricated outputs'))
        writer.writerow(dict(stage='G2',parent='30 new parents NOT_CREATED',noise='2 each NOT_CREATED',method='all',seed='not opened',status='NOT_RUN_PREREQUISITE',reason='G0/G1 not completed; no G2 grant'))
    compiled=[]
    for folder in ['src','tools','tests']:
        for path in (ROOT/folder).rglob('*.py'):
            ast.parse(path.read_text());compiled.append(str(path.relative_to(ROOT)))
    tests=[]
    for path in (ROOT/'results/jobs').glob('*test*/RECEIPT.json'):
        rec=json.loads(path.read_text())
        if rec['exit_code']==0:
            match=re.search(r'Ran (\d+) tests?',(path.parent/'OUTPUT.log').read_text())
            if match:tests.append((path.stat().st_mtime,path,int(match[1])))
    latest=sorted(tests)[-1]
    dump('CODE_CHECK.json',dict(status='syntax_pass',sources=compiled,no_pde=True,
        test_receipt=str(latest[1].relative_to(ROOT)),unittest_count=latest[2]))
    print('Actual mesh figure, costs, NOT_RUN rows and stage decisions saved. No PDE.')

if __name__=='__main__':main()
