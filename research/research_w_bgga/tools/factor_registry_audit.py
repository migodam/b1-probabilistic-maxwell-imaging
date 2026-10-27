"""Check factor replacement on actual exposed observations, with no new PDE.

This is a numerical/integration regression after all online methods were sealed.
It does not change any method, assimilate the outer held-out channels, or add a
new acquisition. The prior is the explicit old finite counting measure.
"""
from pathlib import Path
import argparse,json
import numpy as np
from bgga.intervals import FinitePriorIntervals, LogLikelihoodBounds
from g1_replay import audited_inputs, _load_keys

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--old-bank-root',type=Path,default=None)
    args=parser.parse_args()
    seal = ROOT/'results/g1_replay_001/ALL_ONLINE_OUTPUTS_SEALED.json'
    if not seal.is_file():
        raise RuntimeError('all online outputs must precede this audit')
    inputs = audited_inputs(json.loads((ROOT/'research/delegated/source_audit/CACHE_MAP.json').read_text()),args.old_bank_root)
    payloads = {key: _load_keys(v['path'], v['keys']) for key,v in inputs.items()}
    prior = payloads['FINITE_PRIOR_GRID']['prior_weight']
    bank = payloads['FINITE_FULL_RESPONSE_CACHE']['data'][:, :2, :]
    observations = payloads['FINITE_24_DEVELOPMENT_CASES']['data'][:, :, :2, :]
    rows = []
    for p in range(12):
        for n in range(2):
            log_l = -0.5*np.sum(np.abs((bank-observations[p,n])/0.01)**2, axis=2)
            seq = FinitePriorIntervals(prior)
            ids = [f'bgga.parent{p:02d}.noise{n}.tx{k}.rx0-15.f1' for k in range(2)]
            seq.set_factor(ids[0], LogLikelihoodBounds.exact(log_l[:,0]))
            # Refining/recomputing an existing observation replaces its version.
            first = seq.posterior().joint_weights.copy()
            seq.set_factor(ids[0], LogLikelihoodBounds.exact(log_l[:,0]), version=2)
            refresh_error = float(np.max(np.abs(first-seq.posterior().joint_weights)))
            seq.set_factor(ids[1], LogLikelihoodBounds.exact(log_l[:,1]))
            batch = FinitePriorIntervals(prior)
            batch.set_factor(f'batch.p{p}.n{n}', LogLikelihoodBounds.exact(log_l.sum(axis=1)))
            err = float(np.max(np.abs(seq.posterior().joint_weights-batch.posterior().joint_weights)))
            rejected = False
            try:
                seq.set_factor('illegal_duplicate', LogLikelihoodBounds.exact(log_l[:,0]), observation_id=ids[0])
            except ValueError:
                rejected = True
            assert refresh_error == 0 and err < 1e-12 and rejected
            rows.append(dict(parent=p, noise=n, observation_factor_ids=ids,
                versions=[2,1], refresh_changes_posterior=refresh_error,
                batch_incremental_max_difference=err, duplicate_observation_rejected=rejected))
    result = dict(scope='STATIC_BGGA_V1', status='ACTUAL_CACHED_OBSERVATION_FACTOR_AUDIT_PASS',
        rows=rows, distinct_observation_ids=48, data_dimensions_per_factor=32,
        data_used='same already exposed TX0/TX1 observations, partitioned once',
        no_new_acquisition=True, outer_TX2_3='OFFLINE_PREDICTIVE_SCORING_ONLY_NEVER_ASSIMILATED',
        M4_refinement='versions replace one observation-model factor; queries are not new data',
        acquisition_panel='NOT_RUN', new_full_FE=0, new_complex_RHS=0,
        interpretation='finite numerical factor bookkeeping; not a new Bayesian or Maxwell theorem')
    (ROOT/'FACTOR_HISTORY.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(status=result['status'], cases=len(rows),
        max_batch_error=max(r['batch_incremental_max_difference'] for r in rows))))

if __name__=='__main__': main()
