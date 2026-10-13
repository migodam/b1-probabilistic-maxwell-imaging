"""Read-only derived checks and diagnostics; zero full physics evaluations."""
from pathlib import Path
import json,time
import numpy as np
from scipy.special import ndtr,logsumexp
from src.b1.inference import log_likelihood
from src.b1.budget import Ledger
from src.b1.pilot import dump
start=time.perf_counter();root=Path.cwd();c=json.loads(Path('configs/protocol.json').read_text());led=Ledger('results',384);b=dict(np.load('results/banks/order7.npz'));models=np.load('results/LINEAR_MODELS.npz');records=[]
prior_scale=np.array([.6,.12])/np.sqrt(12)
for name in ['distinguishable','competitive','information_limited']:
 obs=np.load(f'results/cases/{name}/posterior.npz');sigma=float(obs['noise_real_sd']);spectra=[]
 for h in range(3):
  K=models['J_full'][h]*prior_scale/sigma;sv=np.linalg.svd(K,compute_uv=False);spectra.append(dict(H=h,prior_scale_normalized_likelihood_curvature=(sv*sv).tolist(),interpretation='local central-FD metric; uniform-prior scale, not Gaussian posterior covariance'))
 ref=[]
 for h in range(3):
  ids=np.flatnonzero(b['H']==h);i=ids[np.argmin(np.sum((b['theta'][ids]-[2.5,.1])**2,axis=1))];ref.append(b['F'][i])
 pairs=[]
 for i,j in [(0,1),(0,2),(1,2)]:
  distance=np.linalg.norm(ref[i]-ref[j])/sigma;pairs.append(dict(H0=i,H1=j,whitened_full_secant_distance=float(distance),equal_prior_fixed_material_Bayes_error=float(ndtr(-distance/2)),scope='fixed permittivity 2.5+0.1i pairs; NOT marginalized topology error'))
 records.append(dict(case=name,metric=spectra,finite_pairs=pairs))
dump('results/MATERIAL_OBSERVABILITY_AND_SECANTS.json',dict(records=records,new_physics=0,already_paid_derivative_calls=12))
# Negative control: omit the candidate-dependent determinant of gain covariance.
ydata=np.load('results/cases/gain_marginalized/posterior.npz');y=ydata['data'];sigma=float(ydata['noise_real_sd']);L,mah=log_likelihood(b['F'],y,sigma,.15)
ff=np.sum(b['F']**2,axis=1);badL=L+.5*np.log1p(.15**2*ff/sigma**2);pi=ydata['topology_prior']
def event_mass(loglik):
 z=np.array([logsumexp(loglik[b['H']==h]+np.log(b['weights'][b['H']==h])) for h in range(3)]);return np.exp(np.log(pi)+z-logsumexp(np.log(pi)+z))
dump('results/GAIN_LOGDET_NEGATIVE_CONTROL.json',dict(correct_gain_covariance_mass_order7=event_mass(L),omitted_logdet_WRONG_mass_order7=event_mass(badL),probability_difference_max=float(max(abs(event_mass(L)-event_mass(badL)))),new_physics=0,note='wrong target; no standalone accuracy claim, order7 original reference unresolved for this case'))
# Image/remaining mass and all reservations independently checked.
rows=led.rows();begin=[r for r in rows if r['event']=='BEGIN'];end=[r for r in rows if r['event']=='END'];flags=[]
for p in sorted(Path('results/cases').glob('*/probabilities.json')):
 d=json.loads(p.read_text());probs=np.asarray(d['result']['topology_probabilities']);k=d['top_k'];total=sum(q['posterior_mass'] for q in k['items'])+k['remaining_modeled_mass'];flags.append(dict(case=d['case'],probability_sum=float(probs.sum()),topk_with_remaining_sum=float(total),pass_sum=abs(total-1)<1e-10 and abs(probs.sum()-1)<1e-10))
summary=dict(status='PASS' if all(f['pass_sum'] for f in flags) and len(begin)==len(end)==384 else 'FAIL',cases=flags,physics=led.summary(),all_ids_complete=[r['id'] for r in begin]==list(range(1,385)) and sorted(r['id'] for r in end)==list(range(1,385)),wall_seconds=time.perf_counter()-start,new_physics=0)
dump('results/FINAL_DERIVED_AUDIT.json',summary);print(json.dumps(dict(status=summary['status'],cases=len(flags),new_physics=0,wall_seconds=summary['wall_seconds'])))
