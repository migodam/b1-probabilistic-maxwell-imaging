"""Analysis-only summaries and figures; all probabilities carry reference qualification."""
from pathlib import Path
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=Path(__file__).resolve().parents[1]
def load(p):return json.loads(p.read_text())
def main():
    rows={}
    for folder in ['v1_failed_preserved/physical','v3_strata','remaining_spaces']:
        for p in (OUT/folder).glob('*_RESULT.json'):
            d=load(p);d['source']=str(p.relative_to(OUT));rows[d['observation']]=d
    table=[]
    for si in range(3):
        for pi in range(2):
            for ni,mult in enumerate([1,.5,.25]):
                key=f'S{si}_P{pi}_N{ni}';d=rows.get(key,{});co=d.get('original_convergence',{});cr=d.get('cross_coordinate_agreement',{})
                table.append(dict(observation=key,space=si,parent=pi,noise_multiplier=mult,status=d.get('status',d.get('qualification','NOT_COMPLETED')),original_Z_change=co.get('relative_Z_change'),original_material_mean_change_SD=co.get('material_mean_error_over_SD'),response_Z_difference=cr.get('relative_Z_change'),response_material_mean_difference_SD=cr.get('material_mean_error_over_SD'),RTO_successes=d.get('RTO_successes'),source=d.get('source'),scale_up='HOLD'))
    # Append later rank2 cross-coordinate evidence without replacing original reference diagnostics.
    for row in table:
        p=OUT/f"rank2_higher_field/{row['observation']}_RESULT.json"
        if p.exists():
            d=load(p);co=d['cross_coordinate_agreement'];row['response_Z_difference']=co['relative_Z_change'];row['response_material_mean_difference_SD']=co['material_mean_error_over_SD']
        row['numerical_reference_moment_gate']=bool(row['original_Z_change'] is not None and row['original_Z_change']<=.01 and row['original_material_mean_change_SD']<=.01)
        row['coverage_gate']='NOT_CERTIFIED'
    (OUT/'OBSERVATION_SUMMARY.json').write_text(json.dumps(table,indent=2)+'\n')
    with (OUT/'OBSERVATION_SUMMARY.csv').open('w') as f:w=csv.DictWriter(f,list(table[0]));w.writeheader();w.writerows(table)
    fig,ax=plt.subplots(1,3,figsize=(13,4))
    direct=load(OUT/'v3_strata/HOST_DIRECT_READBACK.json');xx=np.array([1,.5,.25]);mass=[d['host_mass'] for d in direct]
    ax[0].plot(xx,mass,'o-');ax[0].axhline(.2208604823,color='gray',ls='--',label='prior host mass');ax[0].set(xscale='log',xlabel='noise / base sigma',ylabel='host posterior mass',title='Positive-mass rank0 stratum');ax[0].legend(fontsize=8)
    for pi,color in [(0,'tab:orange'),(1,'tab:blue')]:
        ts=[r for r in table if r['space']==1 and r['parent']==pi]
        ax[1].plot([r['noise_multiplier'] for r in ts],[r['original_material_mean_change_SD'] for r in ts],'o-',color=color,label=f'parent {pi}: mode partition')
        ax[2].plot([r['noise_multiplier'] for r in ts],[r['response_material_mean_difference_SD'] for r in ts],'o-',color=color,label=f'parent {pi}')
    ax[1].plot(xx,[d['convergence']['material_mean_error_over_SD'] for d in direct],'s--',color='tab:red',label='parent0: direct stratum readback')
    for a in ax[1:]:a.axhline(.01,color='black',ls=':');a.set(xscale='log',yscale='log',xlabel='noise / base sigma',ylabel='material mean error / posterior SD');a.legend(fontsize=7)
    ax[1].set_title('Original-coordinate convergence');ax[2].set_title('Response vs original (qualified limits)');fig.tight_layout();fig.savefig(OUT/'figures/A_PHYSICAL_CONVERGENCE.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(10,4))
    diagnostics=[]
    for ni,mult in enumerate([1,.5,.25]):
        p=OUT/f'v3_strata/S1_P1_N{ni}_RESPONSE_9.json'
        if not p.exists():continue
        d=load(p);rec=d['diagnostics']['records'];l=np.array([r['log_mass'] for r in rec]);w=np.exp(l-l.max());w/=w.sum();a=np.array([r['a'] for r in rec]);s=np.array([r['s'] for r in rec]);order=np.argsort(a);mu=w@a;sd=np.sqrt(w@((a-mu)**2))
        ax[0].step((a[order]-mu)/sd,np.cumsum(w[order]),where='post',label=f'noise={mult}')
        ax[1].scatter(a,s,s=7,c=w,cmap='viridis',alpha=.5)
        diagnostics.append(dict(observation=f'S1_P1_N{ni}',weak_mean=mu,weak_sd=sd,selected_response_mean=float(w@s),selected_response_variance=float(w@((s-w@s)**2)),mean_normal_residual2=float(w@np.array([r['normal_residual2'] for r in rec])),inverse_failures=d['diagnostics']['inverse_failures'],observed_weight_range_is_not_certificate=True))
    ax[0].set(title='Regular transition weak marginal',xlabel='standardized field coordinate',ylabel='numerical CDF');ax[0].legend();ax[1].set(title='Actual weak/selected-response nodes',xlabel='field coordinate',ylabel='selected normalized residual s');fig.tight_layout();fig.savefig(OUT/'figures/A_WEAK_RESPONSE_COORDINATES.png',dpi=170);plt.close(fig)
    (OUT/'RESPONSE_DIAGNOSTICS.json').write_text(json.dumps(diagnostics,indent=2)+'\n')
    print(json.dumps(dict(logical_observations=len(table),with_result=sum(r['source'] is not None for r in table),sources=len(rows))))
if __name__=='__main__':main()
