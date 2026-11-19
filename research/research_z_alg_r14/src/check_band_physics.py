"""Independent non-grid forwards for the saved known-gain continuum witness."""

import json
import os
from time import perf_counter, process_time

for name in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[name] = "1"

from flint import arb
import numpy as np
from scipy.special import j0
from scipy.stats import norm

from continuous_band import Physical, pt
from data14 import ROOT, FrequencyTM, load_cases
from freq_model import Model

WITNESS = ROOT / "results/continuous_band/fresh6-16_known_coherent.json"
DERIVATIVE_POINTS = (-.0197316, -.01461879, -.0034129, .00678137, .0173311, .0196252)


def midpoint(matrix):
    return np.array([[complex(float(matrix[i,j].real.mid()), float(matrix[i,j].imag.mid()))
                      for j in range(matrix.ncols())] for i in range(matrix.nrows())])


def relative(actual, reference):
    return float(np.linalg.norm(actual-reference)/np.linalg.norm(reference))


class Independent:
    def __init__(self, geometry, materials):
        self.geometry = geometry
        self.k0 = geometry["k"]
        self.sources = np.array(geometry["sources"])
        self.materials = materials
        self.cache = {}
        self.fact = self.rhs = 0
        phi = 2*np.pi*np.arange(256)/256
        directions = np.column_stack((np.cos(phi), np.sin(phi)))
        self.phase = directions@self.sources.T
        self.source_distances = np.linalg.norm(self.sources[:, None]-self.sources[None, :],axis=-1)

    def response(self, offset):
        if offset not in self.cache:
            k = self.k0*(1+offset)
            model = Model(**dict(self.geometry,k=k))
            # Keep the raw R2 model and normalize incident/result amplitude here.
            rows = np.array([model.forward(m,independent=True)["y"][0] for m in self.materials])
            self.fact += model.counts["LU"]
            self.rhs += model.counts["actual_complex_RHS"]
            self.cache[offset] = rows*np.sqrt(k/self.k0)
        return self.cache[offset]

    def gram(self, offset):
        # Periodic angular quadrature of source far-field patterns.
        k = self.k0*(1+offset)
        phase = np.exp(1j*k*self.phase)
        power = phase.conj().T@phase/len(phase)
        return power, j0(k*self.source_distances)


def inspect_point(independent, offset, terminal, threshold):
    rows = independent.response(offset)
    difference = rows[0]-rows[1]
    gram, reference_gram = independent.gram(offset)
    direction = np.linalg.solve(gram,difference.conj())
    optimized_square = float(np.real(difference@direction))
    coefficients = direction*np.sqrt(3/optimized_square)
    fixed = float(np.linalg.norm(difference))
    coherent = float(abs(difference@coefficients))
    return dict(frequency_offset=offset, terminal_center=terminal[0], terminal_half_width=terminal[1],
                fixed=fixed, coherent=coherent, threshold_ratio=coherent/threshold,
                fixed_over_saved_upper=fixed/terminal[2], coherent_over_saved_upper=coherent/terminal[3],
                gram_direct_quadrature_bessel_max_error=float(np.max(abs(gram-reference_gram))),
                gram_smallest_eigenvalue=float(np.linalg.eigvalsh(gram)[0]),
                maximizing_power_error=float(abs(np.vdot(coefficients,gram@coefficients)-3)))


def derivative_checks(independent, case):
    physical = Physical(case)
    r12 = FrequencyTM(case.geometry)
    checks = []
    for offset in DERIVATIVE_POINTS:
        kernels = physical.kernels(pt(offset),arb(0))
        value_derivative = [physical.response(m,kernels,False) for m in independent.materials]
        exact = np.array([midpoint(y)[0] for y,_ in value_derivative])
        first = np.array([midpoint(yp)[0] for _,yp in value_derivative])
        second = np.array([midpoint(physical.response(m,kernels,True))[0]
                           for m in independent.materials])
        legacy = [r12.response(m,offset,derivative=True) for m in independent.materials]
        r12_value = np.array([v[0] for v in legacy])
        r12_first = np.array([v[2] for v in legacy])
        point_runs = []
        for step in (1e-4, 3e-5, 1e-5):
            minus, plus = (independent.response(offset-step), independent.response(offset+step))
            fd = (plus-minus)/(2*step)
            second_fd = (plus-2*independent.response(offset)+minus)/step**2
            point_runs.append(dict(step=step, first_relative=relative(fd,first),
                                   first_max_absolute=float(np.max(abs(fd-first))),
                                   first_pair_relative=relative(fd[0]-fd[1],first[0]-first[1]),
                                   second_relative=relative(second_fd,second),
                                   second_pair_relative=relative(second_fd[0]-second_fd[1],second[0]-second[1])))
        gram, _ = independent.gram(offset)
        checks.append(dict(frequency_offset=offset,
                           r2_vs_arb_relative=relative(independent.response(offset),exact),
                           r12_vs_arb_relative=relative(r12_value,exact),
                           r12_first_vs_arb_relative=relative(r12_first,first),
                           gram_quadrature_vs_arb_max_absolute=float(np.max(abs(gram-midpoint(kernels[3])))),
                           centered_differences=point_runs))
    return checks,dict(r12_factorizations=r12.factorizations, arb_point_inverses=physical.solves,
                       arb_point_rhs=physical.rhs)


def run():
    wall,cpu = perf_counter(),process_time()
    record = json.loads(WITNESS.read_text())
    case = next(c for c in load_cases() if c.case_id==record["case_id"])
    independent = Independent(case.geometry,np.array(record["materials"]))
    threshold = record["threshold_2r"]
    grid = []
    for terminal in record["terminal_intervals"]:
        # Both points lie well inside their saved terminal interval.
        for sign in (-1,1):
            frequency = terminal[0] + sign*terminal[1]/np.sqrt(3)
            grid.append(inspect_point(independent,frequency,terminal,threshold))
    extremes = [inspect_point(independent,sign*record["half_band"],
                              max(record["terminal_intervals"],key=lambda x:sign*x[0]),threshold)
                for sign in (-1,1)]
    derivatives, costs = derivative_checks(independent,case)
    best = max(grid,key=lambda r:r["threshold_ratio"])
    widest_ratio = max(grid,key=lambda r:r["coherent_over_saved_upper"])
    lower_error = float(norm.cdf(-max(t[3] for t in record["terminal_intervals"])/(np.sqrt(2)*case.sigma)))
    result = dict(case=case.case_id, source=str(WITNESS.relative_to(ROOT)),
                  materials=record["materials"], gains=record["gains_zeta_ri"],
                  half_band=record["half_band"], saved_intervals=len(record["terminal_intervals"]),
                  off_grid_count=len(grid), largest_off_grid_gap=best, closest_coherent_envelope=widest_ratio,
                  largest_off_grid_fixed_over_envelope=max(r["fixed_over_saved_upper"] for r in grid),
                  largest_quadrature_bessel_error=max(r["gram_direct_quadrature_bessel_max_error"] for r in grid),
                  smallest_gram_eigenvalue=min(r["gram_smallest_eigenvalue"] for r in grid),
                  largest_maximizing_power_error=max(r["maximizing_power_error"] for r in grid),
                  envelope_failures=sum(max(r["fixed_over_saved_upper"],r["coherent_over_saved_upper"])>1 for r in grid),
                  explicitly_sampled_edges=extremes, derivative_checks=derivatives,
                  gaussian_bayes_error_lower_from_saved_upper=lower_error,
                  source_coefficient_total_power=3.,
                  old_scores_squared_r2=(np.linalg.norm(independent.response(0.)-case.old,axis=1)/case.additive_radius)**2,
                  costs=dict(**costs,r2_factorizations=independent.fact,r2_rhs=independent.rhs,
                             cpu_seconds=process_time()-cpu,wall_seconds=perf_counter()-wall))
    result["old_scores_squared_r2"] = result["old_scores_squared_r2"].tolist()
    (ROOT/"results/band_physics_check.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k not in ("derivative_checks","materials")},indent=2))


if __name__=="__main__":
    run()
