"""Parent final method lock BEFORE either new observation cohort exists.

No truth, observation generation, training, proposal evaluation or physics.
Both independent cohorts commit to exactly the same methods and decisions.
"""
from pathlib import Path
import json,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,Book,write
from p3r4.evaluation import checkpoint_identity

def main():
    for prefix in ('R4L1','R4L2'):
        if next((OUT/'inputs/online').glob(prefix+'_*.npz'),None):raise ValueError('cannot retrospectively freeze observations')
    cfg=OUT/'configs'
    if any((cfg/f'FROZEN_{g}.json').exists() for g in ('LOCKED','CONFIRMATION')):raise FileExistsError('final lock immutable')
    if json.loads((OUT/'models/geometry_weighted_cfm_v1/COMPLETE.json').read_text())['status']!='COMPLETE':raise ValueError('weightedCFM fits must allclose first')
    if json.loads((OUT/'results/physics_transport_paired_v2/CONTROLLER_RESULT.json').read_text())['status']!='COMPLETE':raise ValueError('physical refinement must allclose')
    b=Book('IV','freeze_both_fresh_cohorts',analysis_only=True,max_wall=180)
    status,error='COMPLETE',None
    try:
        methods=[]
        definitions=[('raw_spline','density','density_spline_32768',False),
            ('matched_diag','density','extended_diag_gaussian',False),
            ('tri_maxwell','density','extended_triangular_maxwell',False),
            ('tri_random','density','extended_triangular_random',False),
            ('raw_cfm','density','extended_cfm_matched',True),
            ('anchor_cfm_none','anchor_density','anchor_cfm_none',True),
            ('anchor_cfm_geometry','geometry_weighted_cfm_v1','geometry_cfm_geometry',True),
            ('anchor_cfm_random','geometry_weighted_cfm_v1','geometry_cfm_random',True)]
        for family,folder,prefix,cfm in definitions:
            for seed in (2026102501,2026102502,2026102503):
                job=f'{prefix}_s{seed}';ck=OUT/'models'/folder/job/'BEST.pt';receipt=OUT/'I/jobs'/job/'CLOSED.json'
                if json.loads(receipt.read_text())['status']!='COMPLETE':raise ValueError('closed source fit required: '+job)
                m=dict(id=f'{family}_s{seed}',family=family,seed=seed,checkpoint=str(ck),checkpoint_identity=checkpoint_identity(ck),source_receipt=checkpoint_identity(receipt))
                if cfm:m['sample_nfes']=[32]
                methods.append(m)
        for arm in ('physics','supervised'):
            job=f'transport_v2_raw_spline_{arm}_s2026102501';ck=OUT/'models/physics_transport/paired_v2'/job/'BEST.pt';receipt=OUT/'I/jobs'/job/'CLOSED.json'
            if json.loads(receipt.read_text())['status']!='COMPLETE':raise ValueError('paired diagnostic source required')
            methods.append(dict(id='refined_'+arm,family='refined_'+arm,seed=2026102501,checkpoint=str(ck),checkpoint_identity=checkpoint_identity(ck),source_receipt=checkpoint_identity(receipt)))
        density=dict(status='FROZEN_BEFORE_EITHER_NEW_COHORT',methods=methods,
            method_rule='all declared completed families, no selection by pending DEV CFM image or new test result',K=128,fixed_sample_index=0,sample_seed=20261025031,
            CFM_NFE=32,CFM_NFE_scope='16/32/64 already evaluated on development; fresh inference fixes32 for everyCFM; no exactCFMdensity claim',
            primary_probability=dict(candidate='raw_spline',reference='matched_diag',metric='proper full66 NLL;3seedmeanloss perparent'),
            geometry_comparisons=['tri_maxwell vs tri_random same32kdata','anchor_cfm_geometry vs anchor_cfm_none/random same1024anchorcontexts'],
            weighted_CFM_selection='include completed predefined6models, validation latent-energy favorable but not posterior qualification',
            physics_refinement='all3 replication gates FAILED; only fixedseed1 diagnosticSpline pair retained, no repeatseedrescue',
            probability_law='rawGaussianized dx66; no material-volume factor; allCDF/plateau flags retained',
            full_checks=dict(groups=['raw_spline_s2026102501','anchor_cfm_none_s2026102501__nfe32','anchor_cfm_geometry_s2026102501__nfe32','anchor_cfm_random_s2026102501__nfe32'],
                fixed0='all768independentparents, no residualselection',extended='first16L1parents first64draws for all4groups',max_RHS=56832,
                prospective_amendment='4groups ratherthan initial3 to retain both sameinformation geometrycontrols; unchangedtotalcaps; declaredbeforedata'),
            source_identities=[checkpoint_identity(ROOT/p) for p in ['tools/evaluate_r4_density_frozen_nfe.py','tools/evaluate_r4_density.py','src/p3r4/probability_eval.py','src/p3r4/density.py','src/p3r4/anchor_density.py']],
            scientific_status='conditional models and supervised risk actions; no qualified FULL posterior',new_truth_read=False)
        shared=cfg/'FROZEN_DENSITY_BOTH_COHORTS.json';write(shared,density)
        component=OUT/'AUDIT/risk/fresh_freeze_v3/RISK_COMPONENT_FROZEN.json'
        for cohort,upper in [('locked','LOCKED'),('confirmation','CONFIRMATION')]:
            template=component.parent/f'TOP_LEVEL_{upper}_TEMPLATE.json';top=json.loads(template.read_text())
            top.update(status='FROZEN_BEFORE_NEW_OBSERVATIONS',released=True,parent_VI_component_binding=checkpoint_identity(shared),freeze_time_ns=time.time_ns(),
                source_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),risk_and_probability_targets_separate=True,
                optimization=dict(native_evaluations=96,probit_evaluations=64,initialization='measurement/prior only; noGT',array_codec='stored',fallback='originalfiniteanchor; allfailuresretained'),
                principal_test_intervals='98.3333% CI for3predeclaredprimarycomparisons; othercomparisons exploratory, notglobalFWER',
                primary_geometry_risk='strong_spectral vs strong_random same1024trainingcontexts',additional_risk_comparisons='Spline vsreoptimizedsameclass/Diag/permuted; no assumptionFlownecessary',
                representative='first64L1,fixedSpline_s1+R2raw/clipped+sharedordinarycontrol; same16state/256RHScap;StageIV',adapt_after_L1=False,
                full_probability_expansion='NO further correctedposteriorlargecampaign unless existingprespecifiedESS/weightgate met; finitepilotnotteacher',new_3D=False,new_Mie=False,new_statistics='NOT_RUN; thislock concerns originalP1only')
            dest=cfg/f'FROZEN_{upper}.json';write(dest,top)
            folder=OUT/'results/fresh'/cohort;folder.mkdir(parents=True,exist_ok=True)
            command=[sys.executable,'tools/r4_fresh_risk.py','bind','--cohort',cohort,'--component',str(component),'--top-freeze',str(dest),'--out',str(folder/'BOUND_RISK.json')]
            p=subprocess.run(command,cwd=ROOT)
            if p.returncode:raise RuntimeError('risk binding failed before generation')
            spec=dict(methods=methods,parents=top['parents'],sample_seed=20261025031,
                comparisons=[dict(candidate_family='raw_spline',reference_family='matched_diag',metric='full66_nll'),
                    dict(candidate_family='tri_maxwell',reference_family='tri_random',metric='full66_nll'),
                    dict(candidate_family='anchor_cfm_geometry__nfe32',reference_family='anchor_cfm_random__nfe32',metric='material_energy_score')],
                cohort='NEW_INDEPENDENT_'+upper,preobservations_component=checkpoint_identity(shared))
            write(folder/'DENSITY_METHODS.json',spec)
            p=subprocess.run([sys.executable,'tools/evaluate_r4_density_frozen_nfe.py','freeze','--methods',str(folder/'DENSITY_METHODS.json'),'--out',str(folder/'FROZEN_DENSITY.json')],cwd=ROOT)
            if p.returncode:raise RuntimeError('density binding failed before generation')
        write(OUT/'results/FRESH_METHODS_FROZEN_BEFORE_OBSERVATIONS.json',dict(status='COMPLETE',new_observations=0,new_truth_read=False,new_Maxwell_RHS=0,cohorts=[512,256],same_rule_both=True))
    except BaseException as exc:status,error='FAILED',repr(exc);raise
    finally:b.close(status,error)

if __name__=='__main__':main()
