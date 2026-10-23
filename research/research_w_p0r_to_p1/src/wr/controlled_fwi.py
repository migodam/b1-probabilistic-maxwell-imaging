"""Graded full-wave optimizer controls; known geometry, unknown two properties.

These are actual 2D PDE experiments, not replacements for irregular-wall
online reconstruction. Test labels are public controlled inputs only.
"""
import json
import numpy as np
from wq.meshes import interface_mesh,circle_loop
from .physics import Forward,coefficient_jacobian
from .optimization import trust_region,lbfgsb
from .accounting import ROOT,Book,write,BudgetStop

class ControlModel:
    n=2
    def __init__(self,kind):self.kind=kind
    def shape(self,p):
        if self.kind=='disk':return (np.linalg.norm(p-[0,1.5],axis=1)<.3).astype(float)
        m=(abs(p[:,0])<.7)&(p[:,1]>.6)&(p[:,1]<1.2)
        return m*(.7+.3*np.cos(2*p[:,0])**2 if self.kind=='smooth_block' else 1.)
    def coefficient(self,x,f):return lambda p:1+self.shape(p)*(8*x[0]+1j*.15*x[1]/f)
    def jacobian(self,p,s,x):
        # Two parameter JVPs cost eight RHS, less than twelve receiver solves.
        real=p.jvp(s,lambda pts:8*self.shape(pts)).ravel()
        loss=p.jvp(s,lambda pts:1j*.15*self.shape(pts)/p.freq).ravel()
        return np.c_[real,loss],0.
    def material_metadata(self,x):return dict(model='controlled_known_geometry',shape=self.kind,identity='CONTROLLED_OPTIMIZER_TEST_NOT_UNKNOWN_WALL')

def main():
    book=Book('R');rows=[]
    for kind,truth in [('vacuum',[0.,0.]),('block',[.375,.2]),('disk',[.5,.3]),('smooth_block',[.3,.25])]:
        loops=[circle_loop([0,1.5],.3,256)] if kind=='disk' else [np.array([[-.7,.6],[.7,.6],[.7,1.2],[-.7,1.2]])]
        mesh=interface_mesh(h=.18,loops=loops);model=ControlModel(kind)
        fw=Forward(model,mesh,[1.],book,f'R3/controls/{kind}/truth',phase='R3_control_generation',identity='CONTROLLED_GENERATION')
        y,_=fw.evaluate(np.asarray(truth))
        for method in ['TR','LBFGSB']:
            ff=Forward(model,mesh,[1.],book,f'R3/controls/{kind}/{method}',phase='R3_control_inversion',identity='MEASUREMENT_ONLY_WITH_KNOWN_CONTROL_GEOMETRY',retain_current=True)
            out=ff.folder;hist=[]
            def cb(row):hist.append(row);write(out/'ITERATIONS.json',hist)
            x0=np.array([.1,.1]);solver=trust_region if method=='TR' else lbfgsb
            kw=dict(max_iter=14,max_calls=32,callback=cb,fatal_exceptions=(BudgetStop,))
            if method=='TR':kw.update(radius=.15,l2_radius=False,gtol=1e-6)
            res=solver(ff.residual(y),x0,np.zeros(2),np.ones(2),**kw)
            d,_=ff.evaluate(res.x)
            row=dict(kind=kind,method=method,d_mean=float(np.linalg.norm(d-y)/.01),projected_gradient=float(res.projected_gradient),success=bool(res.success),stop=str(res.message),iterations=res.nit,recovered_properties=res.x,identity='CONTROLLED_OPTIMIZER_TEST_NOT_UNKNOWN_WALL')
            rows.append(row);write(ROOT/'results/R3/controls/RESULTS.json',rows);print(json.dumps(row,default=lambda v:v.tolist()),flush=True)
    write(ROOT/'results/R3/controls/COMPLETE.json',dict(totals=book.totals()))
if __name__=='__main__':main()
