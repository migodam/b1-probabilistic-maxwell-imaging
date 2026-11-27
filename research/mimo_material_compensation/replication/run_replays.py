"""Run representative evidence cold, with only this study and archived dependencies."""
from argparse import ArgumentParser
from contextlib import nullcontext
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter, process_time
import importlib.metadata
import json
import os
import shutil
import subprocess
import sys

from restore_prior import restore

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent
BRANCH=Path("research/mimo_material_compensation")
LOCAL_DATA={
    "paths":("n32_c01_paths_ball.json","n32_c01_paths_ball.npz","n32_c01_reduction.npz"),
    "environment":("wave_96_2.0_1.npz","cover_96_2.0_wave_96_2.0_1.json",
                   "verified_responses_96_2.0_wave_96_2.0_1_400.npz",
                   "verified_responses_96_2.0_wave_96_2.0_1_460.npz"),
    "budgets":("freeze.json",),
}
RERUN_PATH=r"""
import json
from pathlib import Path
from threadpoolctl import threadpool_limits
import replay
with threadpool_limits(1):
    rows=replay.inspect(replay.HERE/'n32_c01_paths_ball.json')
    assert rows and all(r['residual_upper_over_r']<1 and r['segments']>0 for r in rows)
    print(json.dumps(rows))
"""
RERUN_TOPOLOGY=r"""
import json
from pathlib import Path
import inspect
import dual_cie
from verify import *
f=HERE/'loop'/'fresh_bg_b'
d=load(f/'design_candidates.npz',('roi',))
roi=d['roi']
engine=PortEngine()
old=load(f/'initial_A.npz')
new=load(f/'new_observed.npz')
spec=json.loads((f/'design_frozen.json').read_text())
rows=[]
with threadpool_limits(1):
    w=VerifiedWave(dict(engine.op,E0=engine.b,S=engine.c))
    for layout in ('B_interleaved','compensation'):
        a=load(f/f'budget_replay_{layout}.npz')
        cfg=spec['actions'][layout]
        center=a['center']
        assert np.all(center[roi]==0.) and center.min()>=0 and center.max()<=4.5
        for key in ('center','fixed_constant_ROI_0.500001','shape_free_ROI_mean_0.500001'):
            x=a[key]
            st=engine.state(x)
            response,error,proof=audit_response(w,x,st['j'])
            old_row=audit_projection(response,error,A,np.arange(8,16),old['Y'],float(old['radius']))
            new_row=audit_projection(response,error,cfg['tx'],cfg['rx'],new[layout+'_Y'],float(new[layout+'_r']))
            exterior=up((fro(x[~roi]-center[~roi])+rd(2)*(fro(x[~roi])+fro(center[~roi])))*(2/32))
            mean_lb=float(np.nextafter(x[roi].mean()*(1-rd(2*roi.sum())),-np.inf))
            assert x.min()>=0 and x.max()<=4.5
            assert old_row['strict_feasible'] and new_row['strict_feasible'] and exterior<.35
            if key!='center':
                assert mean_lb>.5
            if key.startswith('fixed'):
                assert np.all(x[roi]==.500001)
            rows.append(dict(layout=layout,material=key,ROI_mean_lower=mean_lb,exterior_upper=exterior,
                              old_r_upper=old_row['upper_over_r'],new_r_upper=new_row['upper_over_r']))
    provider_paths={k:inspect.getfile(getattr(dual_cie,k)) for k in ('NativeFourierAPCS','Regularizer','_readout')}
print(json.dumps(dict(replayed=rows,CIE_import_paths=provider_paths,physical_calls=dict(engine.counters))))
"""


def run(workspace,artifact):
    start=perf_counter()
    install=restore(workspace)
    study=workspace/BRANCH
    env=dict(os.environ)
    for k in ("PYTHONPATH","PYTHONHOME"):
        env.pop(k,None)
    for k in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","VECLIB_MAXIMUM_THREADS","MKL_NUM_THREADS"):
        env[k]="1"
    env["PYTHONDONTWRITEBYTECODE"]="1"
    env["MPLBACKEND"]="Agg"
    def run_one(name,cmd,mod=None):
        e=dict(env)
        if mod:
            e["PYTHONPATH"]=str(study/mod)
        t,c=perf_counter(),process_time()
        import resource
        before=resource.getrusage(resource.RUSAGE_CHILDREN)
        run=subprocess.run([sys.executable,"-B",*cmd],cwd=workspace,env=e,
                           capture_output=True,text=True)
        after=resource.getrusage(resource.RUSAGE_CHILDREN)
        log=dict(command=["python","-B",*cmd],returncode=run.returncode,
                 wall_s=perf_counter()-t,cpu_s=after.ru_utime-before.ru_utime+after.ru_stime-before.ru_stime,
                 stdout=run.stdout,stderr=run.stderr)
        (artifact/f"{name}.log.json").write_text(json.dumps(log,ensure_ascii=False,indent=2)+"\n")
        run.check_returncode()
        return log

    for sub in ("paths","environment","budgets","topology"):
        dst=study/sub
        dst.mkdir(parents=True)
        for script in (SOURCE/sub).glob("*.py"):
            shutil.copy2(script,dst/script.name)
        for n in LOCAL_DATA.get(sub,()):
            shutil.copy2(SOURCE/sub/n,dst/n)
    for n in ("campaign.json",):
        shutil.copy2(SOURCE/"topology"/n,study/"topology"/n)
    c=Path("topology/loop/fresh_bg_b")
    (study/c).mkdir(parents=True)
    for n in ("initial_A.npz","new_observed.npz","design_candidates.npz","design_frozen.json",
              "budget_replay_B_interleaved.npz","budget_replay_compensation.npz"):
        shutil.copy2(SOURCE/c/n,study/c/n)
    assert not (workspace/"research/z_alg_optical_inertia/blind_uq_mimo_theory/ports/private").exists()
    clocks={}
    logs={}
    def register(name,log):
        logs[name]=log
        clocks[name]={k:log[k] for k in ("wall_s","cpu_s")}
    register("path",run_one("path",["-c",RERUN_PATH],"paths"))
    register("frequency",run_one("frequency",[str(BRANCH/"environment"/"matched_amplitude.py")]))
    register("topology_rank",run_one("topology_rank",[str(BRANCH/"topology"/"rank.py")]))
    register("topology_witnesses",run_one("topology_witnesses",["-c",RERUN_TOPOLOGY],"topology"))
    results=dict(fixed_prior_install=install,workspace=str(workspace),
        packages={pkg:importlib.metadata.version(pkg) for pkg in ("numpy","scipy","python-flint","matplotlib","threadpoolctl")},
        python=sys.version.split()[0],there_are_no_private_old_truth_inputs=True,
        elapsed_wall_s=perf_counter()-start,clocks=clocks,
        path_replay=json.loads(logs["path"]["stdout"]),topology_witnesses=json.loads(logs["topology_witnesses"]["stdout"]))
    frequency=json.loads((study/"environment"/"matched_amplitude.json").read_text())
    ranking=json.loads((study/"topology"/"fixed_pair"/"summary.json").read_text())
    assert frequency["outputs"]["400"]["independent_open_point_pair_upper_over_two_r"]<1
    finite_class=next(row for row in frequency["outputs"]["460"]["budgets"] if row["beta"]==.01)
    assert finite_class["open_lower_over_two_r"]>1 and finite_class["status"]=="all_pairs_excluded"
    assert ranking["cut_count"]==6435 and len(ranking["cases"])==7
    provider_paths=results["topology_witnesses"]["CIE_import_paths"]
    assert all(Path(provider).is_relative_to(workspace) for provider in provider_paths.values())
    (artifact/"native_frequency.json").write_text(json.dumps(frequency,ensure_ascii=False,indent=2)+"\n")
    (artifact/"rank_summary.json").write_text(json.dumps(ranking,ensure_ascii=False,indent=2)+"\n")
    results["topology_cut_count"]=ranking["cut_count"]
    results["topology_object_count"]=len(ranking["cases"])
    (artifact/"replay_receipt.json").write_text(json.dumps(results,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(clocks,ensure_ascii=False,indent=2),flush=True)


if __name__=="__main__":
    p=ArgumentParser()
    p.add_argument("--keep-workspace",type=Path)
    a=p.parse_args()
    artifact=HERE/"replay_receipts"
    artifact.mkdir(exist_ok=True)
    with (nullcontext(str(a.keep_workspace)) if a.keep_workspace else TemporaryDirectory(prefix="mimo-comp-cold-")) as wd:
        workspace=Path(wd).resolve()
        if workspace.exists() and any(workspace.iterdir()):
            raise FileExistsError(f"Expected empty replay workspace: {workspace}")
        workspace.mkdir(parents=True,exist_ok=True)
        run(workspace,artifact)
