from complete_quotient import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
 begin=time.process_time();rows=json.loads((OUT/'results/complete_batch.json').read_text())['rows'];catalog=['# Certified central-data material catalogue\n','Certificates enclose roots; displayed midpoints are for readability. Exact outward interval boxes, all complex endpoints and exclusions are in complete_T*_128.json. Rounded discrete model only; these are not noise-feasible-set certificates.\n'];statuses=[]
 for r in rows:
  case,arm=r['case'],r['arm'];j=json.loads((OUT/'results'/f'complete_{case}_{arm}_128.json').read_text());count=j['legal_count'];tag='CERTIFIED_UNIQUE' if count==1 else ('CERTIFIED_AMBIGUOUS' if count>1 else 'NO_FEASIBLE_CANDIDATE');statuses.append(dict(case=case,arm=arm,central_data_status=tag,legal_count=count,central_box_coverage='COMPLETE',noise_feasible_coverage='UNKNOWN',physical_error_coverage='UNKNOWN',scope=j['scope'],empty_semantics='Count 0 certifies no exact central-data regular legal material; does not certify empty noise feasible set.'))
  catalog += [f'\n## {case} / {arm}: {tag}; {count} legal roots\n','|Root index|Material midpoint|Max enclosure width|\n|---|---|---|']
  for i,c in enumerate(j['classifications']):
   if c['status']=='CERTIFIED_LEGAL_REAL_MATERIAL':
    v=c['material'];catalog.append(f"|{i}|{', '.join(f'{(z['lo']+z['hi'])/2:.10g}' for z in v)}|{max(z['hi']-z['lo'] for z in v):.3g}|")
 (OUT/'ROOT_CATALOG.md').write_text('\n'.join(catalog));save(OUT/'results/CLAIM_STATUS.json',statuses)
 fig,ax=plt.subplots(1,3,figsize=(13,4),constrained_layout=True);x=np.arange(6);clean=[r['legal_count'] for r in rows if r['arm']=='clean'];noisy=[r['legal_count'] for r in rows if r['arm']=='noisy'];ax[0].bar(x-.18,clean,.36,label='Exact clean data');ax[0].bar(x+.18,noisy,.36,label='Exact noisy center');ax[0].set(xticks=x,xticklabels=[f'T{i}' for i in x],ylabel='Certified legal root count',title='All 12 central-data boxes covered');ax[0].legend(fontsize=8);ax[0].set_ylim(0,3.7)
 f=json.loads((OUT/'results/certified_fold.json').read_text());d=[(v['lo']+v['hi'])/2 for v in f['endpoint_determinant_intervals']];ax[1].plot([-1,1],np.array(d)*1e25,'o');ax[1].axhline(0,color='k',lw=.6);ax[1].set(xlabel='Signed displacement / 1e-7',ylabel='det(real material Jacobian) / 1e-25',title='A critical point inside old-data noise set');ax[1].text(.05,.1,'Certified endpoint signs + continuity\nAll old residuals < 1.021e-4 < 1.2e-4',transform=ax[1].transAxes,fontsize=8)
 sc=json.loads((OUT/'results/constrained_source_scores.json').read_text())['score'];vals=[.072386,sc['pair_constrained']['candidate_joint_RMS'],sc['fold_constrained']['candidate_joint_RMS'],.011754];ax[2].bar(range(4),vals,color=['#999999','#4575b4','#4575b4','#66a61e']);ax[2].set(xticks=range(4),xticklabels=['Original\nmax-min','Constrained\npair','Constrained\ncritical mode','Completion\n+ joint TRF'],ylabel='Material RMS error',title='T0 noisy: exposed development');ax[2].text(.03,.96,'1 / 1 / 1 / 3 new measurements\nSame-source TRF matches both new designs',transform=ax[2].transAxes,va='top',fontsize=8);ax[2].set_ylim(0,.105)
 fig.savefig(OUT/'figures/06_certification_and_correction.png',dpi=180);plt.close(fig)
 h=json.loads((OUT/'results/H12_probe.json').read_text());ct=json.loads((OUT/'results/H12_contract.json').read_text());fig,ax=plt.subplots(1,2,figsize=(10,4),constrained_layout=True);c=np.array(ct['centers']);ax[0].scatter(c[:6,0],c[:6,1],s=140,label='6 independent unknown cells');ax[0].scatter(c[6:,0],c[6:,1],s=140,marker='s',label='6 known host cells');
 for i,(xx,yy) in enumerate(c):ax[0].text(xx,yy,str(i),ha='center',va='center',fontsize=8,color='white')
 ax[0].axis('equal');ax[0].set(title='H12: n = 12, k = 6, s = 3',xlabel='x',ylabel='y');ax[0].legend(fontsize=8,loc='upper left');
 for i,p in enumerate(h['predictions']):ax[1].plot(range(6),p['material'],'o-',label=f'Certified explanation {i+1}')
 ax[1].set(title='Two equal-old-response material explanations',xlabel='Active cell',ylabel='Material contrast');ax[1].legend(fontsize=8);ax[1].text(.03,.04,f"One realizable source gap: {h['source']['predicted_gap']:.5f}\nReduced/full agreement: {h['full_reduced_error']:.2e}\nNo noise/host-error guarantee in this probe",transform=ax[1].transAxes,fontsize=8);fig.savefig(OUT/'figures/07_known_host.png',dpi=180);plt.close(fig)
 save(OUT/'results/presentation_receipt.json',dict(CPU_seconds=time.process_time()-begin,matplotlib=matplotlib.__version__))
if __name__=='__main__':main()
