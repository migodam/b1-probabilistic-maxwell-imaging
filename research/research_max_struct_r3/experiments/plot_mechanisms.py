from pathlib import Path
import json,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
def main():
    wall,cpu=time.perf_counter(),time.process_time()
    d=np.load(ROOT/'results/onecell/arrays.npz');r=json.loads((ROOT/'results/onecell/RECEIPT.json').read_text())
    vec=json.loads((ROOT/'results/vector_scope/RECEIPT.json').read_text());fib=json.loads((ROOT/'results/fiber/RECEIPT.json').read_text())
    fig,axs=plt.subplots(2,2,figsize=(11.8,8.4),layout='constrained');color='#167b82';bad='#b63d32'
    ax=axs[0,0];y=d['y'];target=d['target'].item();ax.plot(y.real,y.imag,color=color,lw=2)
    ax.scatter([target.real],[target.imag],c='#248344',label='true m=1.2 / rational + CIE readout',zorder=4)
    ax.scatter([y[-1].real],[y[-1].imag],c=bad,marker='x',s=60,label='local optimizer: m=20',zorder=4)
    ax.set(xlabel='Re observed scattered field',ylabel='Im observed scattered field',title='One physical disk: a circular response arc');ax.axis('equal');ax.legend(fontsize=8,loc='best');ax.grid(alpha=.15)
    ax=axs[0,1];loss=.5*abs(y-target)**2/abs(target)**2
    ax.plot(d['m'],loss,color=color,lw=2);ax.axvline(r['other_data_stationary_material'],color='#777',ls=':',label='interior maximum 16.8795')
    ax.scatter([20],[r['initial_objective']],c=bad,marker='x',s=60,label='constrained local minimum')
    ax.scatter([1.2],[0],c='#248344',label='global data minimum')
    ax.set(xlabel='real material contrast m',ylabel='normalized data objective',title='Radiation passivity does not remove a box trap');ax.legend(fontsize=8);ax.grid(alpha=.15)
    ax=axs[1,0]
    for family,sty,col in [('dyad','o-',color),('isotropic','s-',bad)]:
        vals=vec['response_change_spectra'][family]['relative_to_first'];ax.semilogy([1,2,3],vals,sty,color=col,label=family)
    ax.set(xticks=[1,2,3],xlabel='response-change singular value index',ylabel='relative singular value',title='3D: rank-one dyad vs rank-three isotropic cell');ax.legend();ax.grid(alpha=.15)
    ax.text(.02,.03,'8-cell finite vector DDA; no imaging/continuum claim',transform=ax.transAxes,fontsize=8)
    ax=axs[1,1];ax.bar(['old data','common added loss'],[fib['old_pair_gap'],fib['new_pair_gap']],color=[color,bad]);ax.set_yscale('log')
    ax.set(ylabel='absolute full-port pair gap',title='Existing 6-cell material fiber: altered data separate')
    ax.text(.03,.05,'Old near-equality is numerical, not a new exact certificate',transform=ax.transAxes,fontsize=8);ax.grid(axis='y',alpha=.15)
    fig.suptitle('MAX-STRUCT-R3 | formula-derived physical mechanism controls',fontsize=15)
    out=ROOT/'figures';out.mkdir(exist_ok=True)
    fig.savefig(out/'PHYSICAL_MECHANISMS.png',dpi=170);fig.savefig(out/'PHYSICAL_MECHANISMS.pdf');plt.close(fig)
    rec=dict(status='COMPLETED',input='saved onecell/vector_scope/fiber arrays and receipts',physical_LU=0,physical_RHS=0,
        wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,scope='plots only, no solver imports')
    (ROOT/'audit/MECHANISM_PLOT_RECEIPT.json').write_text(json.dumps(rec,indent=2))
if __name__=='__main__':main()
