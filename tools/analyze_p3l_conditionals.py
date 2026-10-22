"""Zero-new-physics analysis of actual, already paid conditional weights."""
from pathlib import Path
import sys,json,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,write,save

def main():
    b=Book('B','conditional_GN_departure_audit',max_wall=600);status='COMPLETE';error=None
    try:
        report=json.loads((OUT/'B/CONDITIONAL_PILOT_RESULTS.json').read_text());rows=[]
        for row in report['rows']:
            folder=OUT/f'B/conditionals/{row["parent"]}_d{row["weak_dimension"]}_{row["context"]}'
            with np.load(folder/'CONDITIONAL_REFERENCE.npz') as z:mean=z['mean_b'];H=z['precision_b']
            with np.load(folder/'IMPORTANCE.npz') as z:delta=z['strong_samples']-mean;lw=z['logweights']
            ginf=row['strong_gradient_inf'];lambda_min=float(np.linalg.eigvalsh(H)[0]);s=len(mean)
            # Completing the square in the local GN quadratic adds .5*g'C*g.
            # Without the saved vector g its rigorous norm upper bound suffices.
            centering_logZ_upper=.5*s*ginf**2/max(lambda_min,np.finfo(float).tiny)
            departure=lw-row['GN_logZ']
            linear_abs_bound=ginf*np.abs(delta).sum(axis=1)
            non_GN_abs_lower=np.maximum(0.,np.abs(departure)-linear_abs_bound)
            extra=dict(quadratic_centering_logZ_correction_upper_bound=centering_logZ_upper,
                       GN_logweight_departure_median=float(np.median(np.abs(departure))),
                       beyond_local_linear_gradient_absolute_departure_lower_median=float(np.median(non_GN_abs_lower)),
                       beyond_linear_departure_min=float(np.min(non_GN_abs_lower)),
                       maximum_possible_linear_gradient_effect=float(np.max(linear_abs_bound)),
                       interpretation='Departures from frozen GN quadratic include non-GN second curvature and higher terms; cannot assign all to Maxwell feedback',
                       original_GN_normalizer_is_not_stationary_Laplace=row['mode_converged'] is False)
            save(folder/'GN_DEPARTURE_AUDIT.npz',logweight_minus_GN_logZ=departure,linear_gradient_abs_upper=linear_abs_bound,beyond_local_linear_abs_lower=non_GN_abs_lower)
            rows.append({k:row[k] for k in ['parent','weak_dimension','context','importance_ESS','mode_converged']}|extra)
        write(OUT/'B/GN_DEPARTURE_AUDIT.json',dict(new_physics_RHS=0,old_results_unmodified=True,not_a_new_posterior_reference=True,rows=rows,
            conclusion='Small gradient alone does not fix conditional importance collapse. Actual full target is strongly nonquadratic over the declared GN Gaussian support. Gaussian decoder means remain reference quantities, not validated Bayes images.'))
        print(json.dumps(dict(rows=len(rows),centering_upper_below_01=sum(r['quadratic_centering_logZ_correction_upper_bound']<.01 for r in rows),median_beyond_linear_departure=float(np.median([r['beyond_local_linear_gradient_absolute_departure_lower_median'] for r in rows])))))
    except BaseException as e:status='ERROR';error=repr(e);raise
    finally:b.close(status,error)

if __name__=='__main__':main()
