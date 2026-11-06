from pathlib import Path
import csv,collections,json
import numpy as np
P=Path(__file__).resolve().parent
def rows(name):return list(csv.DictReader((P/name).open()))
s=rows('R8_SELECTION_AND_CORRECTOR_630.csv');w=rows('R8_INITIALIZER_SELECTION_SURVIVAL.csv');summary=[]
for init in ['uniform','nn']:
 for arm in ['cie_single','cie_dual','global_dual','a_material','a_plus','a_plus_sourcewise','physical_lbfgs']:
  r=[x for x in s if x['initialization']==init and x['arm']==arm and float(x['amplitude'])>=12]
  if not r:continue
  q=[x for x in w if x['initialization']==init and x['arm']==arm and float(x['amplitude'])>=12]
  z=dict(initialization=init,arm=arm,n=len(r),selected_stages=dict(collections.Counter(x['selected_auxiliary_stage'] for x in r)))
  z['Phi_improved_by_corrector']=sum(float(x['corrector_Phi_reduction'])>0 for x in r)
  z['numerical_success']=sum(x['corrector_success']=='True' for x in r)
  z['numerical_success_material_fail']=sum(x['corrector_success']=='True' and x['material_image_success']=='False' for x in r)
  z['material_worse_than_initial']=sum(float(x['material_nrmse'])>float(x['initial_material_nrmse']) for x in r)
  z['lower_prePhi_but_both_finalPhi_and_material_worse_NNdirect']=sum(x['lower_pre_Phi_but_final_Phi_and_material_worse_NNdirect']=='True' for x in q)
  z['median_pre_to_final_Phi_ratio']=float(np.median([float(x['selected_pre_Phi'])/float(x['original_Phi']) for x in r]))
  summary.append(z)
(P/'SELECTION_HIGH_CONTRAST_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
for x in summary:print(json.dumps(x))
allrows=rows('R8_ENDPOINTS_630.csv')
for case in ['r8_test_0000','r8_test_0003','r8_test_0006','r8_test_0027','r8_test_0029']:
 rr=[x for x in allrows if x['case_id']==case and x['initialization']=='uniform']
 print('CASE',case,[(x['arm'],x['amplitude'],x['relative_complex_RMS_noise'],x['material_nrmse'],x['otsu_iou'],x['selected_auxiliary_stage']) for x in rr])
for case in ['r8_test_0027']:
 rr=[x for x in allrows if x['case_id']==case and x['initialization']=='nn']
 print('NNCASE',case,[(x['arm'],x['material_nrmse'],x['original_Phi'],x['initial_original_Phi'],x['selected_auxiliary_stage']) for x in rr])
