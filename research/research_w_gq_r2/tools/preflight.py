from pathlib import Path
import time,json,resource,numpy as np
from r2.accounting import BudgetBook,plain
from r2.physics import Eval,S1_TX
from r2.worlds import micro_wall,target_defs,make_target
ROOT=Path(__file__).resolve().parents[1]
def main():
    out=ROOT/'results/preflight';out.mkdir(exist_ok=False)
    b=BudgetBook();w=micro_wall(0,1/16,0);t=make_target(target_defs()[0][1])
    st=time.process_time();ev=Eval(b,'E1',1,S1_TX,w.x,w.y)
    state=ev.evaluate(w,t,label='finest resolved public mesh benchmark')
    np.savez_compressed(out/'BENCHMARK.npz',data=state.data,u=state.u,mesh_p=ev.p.mesh.p,mesh_t=ev.p.mesh.t,eps=w.epsinf,ell=w.ell,x=w.x,y=w.y)
    r=dict(cpu_seconds=time.process_time()-st,ndofs=ev.p.ndofs,wall_shape=w.shape,factor_bytes=state.timings['factor_bytes'],backward=state.residual,physical_work=b.snapshot(),RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,micro_cells='interfaces exact; not averaged',status='COST_PREFLIGHT_ONLY',target=t.metadata())
    (out/'RESULTS.json').write_text(json.dumps(r,default=plain,indent=2));print(json.dumps(r,default=plain))
if __name__=='__main__':main()
