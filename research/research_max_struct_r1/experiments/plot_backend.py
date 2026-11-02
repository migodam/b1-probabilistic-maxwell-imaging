"""Plot existing sphere controls only; no Maxwell solves or material queries."""
import json,time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]

def main():
    start=time.perf_counter();cpu=time.process_time()
    out=ROOT/'figures';out.mkdir(exist_ok=True)
    receipt=dict(status='RUNNING',LU=0,RHS=0,new_physical_states=0,
                 input='results/backend_qualification/h_0.04_fields.npz,h_0.03_fields.npz')
    try:
        coarse=np.load(ROOT/'results/backend_qualification/h_0.04_fields.npz')
        fine=np.load(ROOT/'results/backend_qualification/h_0.03_fields.npz')
        fig=plt.figure(figsize=(12,3.7),constrained_layout=True)
        ax=fig.add_subplot(131,projection='3d');p=fine['points']
        ax.scatter(*p.T,s=9,c=p[:,2],cmap='viridis',alpha=.8)
        ax.set(xlabel='x / lambda',ylabel='y / lambda',zlabel='z / lambda',title='Fine sphere: 257 dipoles, h=.03 lambda')
        ax.set_box_aspect((1,1,1));ax.view_init(22,35)
        a=fig.add_subplot(132);b=fig.add_subplot(133)
        deg=np.degrees(fine['theta']);ref=np.linalg.norm(fine['mie_far'],axis=1)
        a.plot(deg,ref,'k',lw=2,label='Analytic Mie')
        for data,label in [(coarse,'DDA h=.04'),(fine,'DDA h=.03')]:
            a.plot(deg,np.linalg.norm(data['far'],axis=1),label=label)
            b.plot(deg,np.linalg.norm(data['far']-data['mie_far'],axis=1)/ref.max(),label=label)
        a.set(xlabel='Scattering angle / degrees',ylabel='Far E amplitude / lambda',title='Absolute amplitude, incident E=x, k=z')
        b.set(xlabel='Scattering angle / degrees',ylabel='Vector error / peak Mie amplitude',title='Complex vector error includes phase')
        for ax2 in (a,b):ax2.grid(alpha=.25);ax2.legend(fontsize=8)
        fig.suptitle('Bounded vector DDA sphere control: epsilon=3, radius=.12 lambda; H(curl) NOT_RUN',fontsize=11)
        fig.savefig(out/'BACKEND_SPHERE_CONTROL.png',dpi=170)
        fig.savefig(out/'BACKEND_SPHERE_CONTROL.pdf');plt.close(fig)
        receipt.update(status='COMPLETED',outputs=['BACKEND_SPHERE_CONTROL.png','BACKEND_SPHERE_CONTROL.pdf'])
    except Exception as exc:
        receipt.update(status='FAILED',error=repr(exc));raise
    finally:
        receipt.update(wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu)
        (out/'BACKEND_PLOT_RECEIPT.json').write_text(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
