"""Posthoc material-invariant mask diagnostic. No inference/physics/new main criterion."""
from core import *
import csv
from scipy.ndimage import label
from score_and_plot import topology
rows=[];idx=load_setup()['ROI']
for ds in ['N','P']:
 for i in range(9):
  oid=f'{ds}_{i:03d}';gt=dict(np.load(P/'data/truth'/f'{oid}.npz'))
  for method in ['FULL256','AGENT_A','AGENT_B','AGENT_C']:
   z=dict(np.load(P/'results'/ds/method/oid/'PREDICTION.npz'));gmask=gt['material'].real>=.1+.65;mmask=z['material'].real>=.1+.65;inter=np.sum(gmask&mmask);union=np.sum(gmask|mmask)
   _,holes=topology(mmask.astype(float));_,gth=topology(gmask.astype(float))
   rows.append(dict(dataset=ds,digit=i+1,method=method,real_contrast_threshold=.65,posthoc=True,primary_threshold_unchanged=True,material_mask_IoU=float(inter/max(1,union)),material_mask_Dice=float(2*inter/max(1,gmask.sum()+mmask.sum())),GT_material_holes=gth,reconstructed_material_holes=holes))
with (P/'audit/POSTHOC_MATERIAL_MASK_DIAGNOSTICS.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
write(P/'audit/NATIVE_AMPLITUDE_GAUGE.json',dict(qualification='ANALYTIC_NATIVE_NEW_MODEL_ONLY',identity='a_new=a/c,delta_new=c*delta preserves D=delta diag(a)',conditions='c>0 andbothoccupancy/propertyboundslegal',P1_Full66_global_gauge_claim=False,new_physics_RHS=0,main_metrics_changed=False,secondary_mask='posthoc fixed realmaterialcontrast .65 diagnostic only,notnewpreregisteredprimary'))
