"""Separate affine-decoder and full-Maxwell secants in a completed search.

Only already-persisted Oracle search arrays are read. No solver, truth data,
optimization, or budget object is created.
"""
import json
from pathlib import Path
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main():
    start = time.process_time()
    source = ROOT/'results/B2/task_preserving_2026102202_narrow'
    out = source/'finite_response_diagnostic'
    out.mkdir(exist_ok=True)
    with np.load(source/'BASE.npz') as a:
        jc, jt = a['J'], a['K']
        cal0, task0 = a['cal'], a['task']
        metric = a['material_metric']
    vectors, rows = {}, []
    for name in ('PROPOSAL', 'FINAL'):
        with np.load(source/(name+'.npz')) as a:
            x, cal, task = a['x'], a['cal'], a['task']
        row = dict(endpoint=name, material_RMS=float(np.sqrt(x@metric@x)),
                   box_max_abs_coordinate=float(np.max(abs(x))))
        for label, jac, f, f0 in [('calibration',jc,cal,cal0), ('task8',jt,task,task0)]:
            linear = (jac@x).reshape(f.shape)/.01
            secant = (f-f0)/.01
            remainder = secant-linear
            nl, ns, nr = (float(np.linalg.norm(v)) for v in (linear,secant,remainder))
            row[label] = dict(linear_prediction_norm=nl, full_secant_norm=ns,
                Maxwell_finite_remainder_norm=nr,
                linear_secant_realified_cosine=float(np.real(np.vdot(linear,secant))/(nl*ns)),
                remainder_fraction_of_full_secant=nr/ns,
                identity_error=float(np.linalg.norm(secant-linear-remainder)))
            for key, value in [('linear',linear),('secant',secant),('remainder',remainder)]:
                vectors[f'{name}_{label}_{key}'] = value
        rows.append(row)
    result = dict(identity='ORACLE_SEARCH_POSTHOC_FINITE_RESPONSE_DIAGNOSTIC',
        source=str(source.relative_to(ROOT)), endpoints=rows,
        material_decoder='LinearWorld: epsilon and conductivity affine in the same 108 material coefficients on a fixed support.',
        decoder_secant_remainder='Identically zero for this affine chart, a property of the implemented formula; no extra physics test claimed.',
        inference='At this finite proposal the linear calibration prediction is too optimistic because of Maxwell nonlinearity, not Gaussian latent saturation. This is one finite case, not a universal obstruction.',
        scientific_gate='Final d_cal>1: no calibration witness, notwithstanding d_task>=3.',
        new_rhs=0, source_physics_cost='All endpoints and baseline Jacobians were already paid in B2.',
        pure_analysis_cpu_seconds=time.process_time()-start)
    np.savez_compressed(out/'SECANT_ARRAYS.npz',**vectors)
    (out/'RESULTS.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(result, allow_nan=False))


if __name__ == '__main__':
    main()
