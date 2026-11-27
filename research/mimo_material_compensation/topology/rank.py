"""Exhaustive fixed witness competition. No clean target readings are used."""
from pathlib import Path
import csv

from ports import (A, B, EDGE_I, EDGE_J, INCIDENCE, PRIOR, RX, TX, WATT, PortEngine,
                   cut_id, dump, edges, elapsed, load, now, project, realify, save)
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import rankdata, spearmanr
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
OUT = HERE / "fixed_pair"
CASES = ["n32_c01","n32_c02","n32_c03","n32_c04","n32_c05","n32_c06","n32_c08"]
IA, IB = cut_id(A), cut_id(B)


def cut_power(z):
    return WATT * (INCIDENCE @ abs(edges(z))**2)


def bistatic_geometry(positions, k):
    direction = positions[:16] / np.linalg.norm(positions[:16], axis=1)[:, None]
    q = k*(direction[EDGE_I]+direction[EDGE_J])
    # For a real material, q/-q are equivalent. Favor distinct, widely separated modes.
    separation = np.minimum(
        np.linalg.norm(q[:,None]-q[None,:],axis=2),
        np.linalg.norm(q[:,None]+q[None,:],axis=2))
    similarity = np.exp(-0.5*(separation/(np.pi/2))**2)
    effective_modes = 64**2 / np.einsum("ij,ij->i", INCIDENCE @ similarity, INCIDENCE)
    return effective_modes, q


def task_basis(centers, roi):
    cx, cy = np.mean(centers[roi], axis=0)
    x, y = centers.T
    # Indicator plus two first spatial moments and quadrupole on the registered ROI.
    basis = np.column_stack((np.ones(len(x)),x-cx,y-cy,((x-cx)*(y-cy))))*roi[:,None]
    return np.linalg.qr(basis)[0]


def fisher_scores(j, centers, roi, radius):
    basis = task_basis(centers, roi)
    variables = edges(j)*np.sqrt(WATT)/radius
    old_ids = INCIDENCE[IA].astype(bool)
    old_j = realify(variables[old_ids])
    old_e = old_j[:,~roi]
    prior_curvature = .01*np.mean(np.sum(old_e**2,axis=1))
    def logdet(j_task):
        moment = (j_task.conj()[..., None]*j_task[:, None, :]).real.reshape(len(j_task),-1)
        block = (INCIDENCE @ moment).reshape(-1,4,4)
        return np.linalg.slogdet(np.eye(4)+block)[1]
    j_task = variables@basis
    gram = old_e @ old_e.T
    gram.flat[::len(gram)+1] += max(prior_curvature, 1e-12)
    ridge = old_e.T @ cho_solve(cho_factor(gram,check_finite=False),old_j@basis)
    # Full old-data nuisance-aware predictive Fisher; correction is frozen before new cut.
    return logdet(j_task), logdet(j_task-variables[:,~roi]@ridge)


def exact_conditional_fisher(j, centers, roi, radius):
    """A plus new D-opt information in a four-dimensional ROI linear model."""
    basis = task_basis(centers,roi)
    all_j = realify(edges(j)*np.sqrt(WATT)/radius)
    ri = all_j@basis
    ei = all_j[:,~roi]
    old_edges = np.flatnonzero(INCIDENCE[IA])
    old = np.concatenate((old_edges,old_edges+len(EDGE_I)))
    reg = max(.01*np.mean(np.sum(ei[old]**2,axis=1)),1e-12)
    gram = ei@ei.T
    a_cov = gram[np.ix_(old,old)] + reg*np.eye(len(old))
    factor = cho_factor(a_cov,check_finite=False)
    cross = gram[:,old]
    j_conditioned = ri - cross @ cho_solve(factor,ri[old],check_finite=False)
    cov = (gram-cross@cho_solve(factor,cross.T,check_finite=False))/reg + np.eye(len(ri))
    static = np.eye(4) + reg*ri[old].T@cho_solve(factor,ri[old],check_finite=False)
    values = np.empty(len(TX))
    for cut,incidence in enumerate(INCIDENCE):
        used = np.flatnonzero(incidence)
        ix = np.concatenate((used,used+len(EDGE_I)))
        local = j_conditioned[ix]
        solved = cho_solve(cho_factor(cov[np.ix_(ix,ix)],check_finite=False),local,check_finite=False)
        values[cut] = np.linalg.slogdet(static+local.T@solved)[1]
    return values


def fairness(engine, material, delta, radius, rng):
    full = engine.full_response(material)
    bare_delta,bare_total,bare_cost = engine.independent_response(material)
    b_ind = np.ix_(np.setdiff1d(np.arange(16),B), B)
    a_test = project(full["delta"], A)
    ar_test = project(full["delta"], np.arange(8,16))
    direct = engine.state(material, columns=A)["delta"][8:16] * np.sqrt(WATT)
    a_edges = set(zip(EDGE_I[INCIDENCE[IA]>0], EDGE_J[INCIDENCE[IA]>0]))
    b_edges = set(zip(EDGE_I[INCIDENCE[IB]>0], EDGE_J[INCIDENCE[IB]>0]))
    c_clean = delta[:16,:16]/4
    b_clean = delta[b_ind]*np.sqrt(WATT)
    c_sym = (c_clean[b_ind]+c_clean.T[b_ind])/np.sqrt(2)
    n=4096
    rc,tb = b_ind
    # Independent Gaussian channel noise. Reciprocal averaging consumes both reads.
    noise_1=(rng.normal(size=(n,8,8))+1j*rng.normal(size=(n,8,8)))/np.sqrt(2)
    noise_2=(rng.normal(size=(n,8,8))+1j*rng.normal(size=(n,8,8)))/np.sqrt(2)
    noise_av=(noise_1+noise_2)/np.sqrt(2)
    empirical_cov=np.cov(noise_av.reshape(n,-1).T)
    rel = lambda a,b: float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-300))
    counts=dict(A=dict(complex_reads=64,reciprocal_complex=64,real_dof_cap=128,pulses=8),
                B=dict(complex_reads=64,reciprocal_complex=64,real_dof_cap=128,pulses=8),
                A_and_B=dict(complex_reads=128,reciprocal_complex=len(a_edges|b_edges),real_dof_cap=2*len(a_edges|b_edges),pulses=16),
                C=dict(complex_reads=256,reciprocal_complex=136,real_dof_cap=272,pulses=16))
    z=dict(physical_ports=24, switchable=list(range(16)), permanently_terminated=list(range(16,24)),
           loads_ohm=engine.zref, frequency_MHz=400,grid=[32,32],counts=counts,
           A_B_reciprocal_overlap=len(a_edges&b_edges),
           total_transpose_reciprocity_rel=rel(full["total"],full["total"].T),
           calibrated_transpose_reciprocity_rel=rel(full["delta"],full["delta"].T),
           A_role_reversal_rel=rel(a_test,ar_test.T),A_direct_vs_projection_rel=rel(direct,a_test),
           bare_VIE_loaded_total_rel=rel(bare_total,full["total"]),
           bare_VIE_loaded_calibrated_rel=rel(bare_delta,full["delta"]),
           original_released_A_coefficients_rel=[rel(engine.b[:,:8]*np.sqrt(WATT),engine.op["E0"]),
                                                 rel(engine.c[8:16],engine.op["S"])],
           unit_power_state_relative_defect=full["state_defect"],
           gaussian_C_symmetrized_signal_equals_B_rel=rel(c_sym,b_clean),
           gaussian_C_channel_variance=[float(np.mean(abs(noise_1)**2)),float(np.mean(abs(noise_2)**2))],
           gaussian_C_to_B_empirical_variance=float(np.trace(empirical_cov).real/64),
           gaussian_C_to_B_empirical_max_offdiagonal=float(np.max(abs(empirical_cov-np.diag(np.diag(empirical_cov))))),
           reference_axis="full S a->b, S=I-2*Zref*(Z_A+Zref*I)^-1, delta=S-S_empty",
           frobenius_radius_A_B=radius,frobenius_radius_C=2*radius,
           source_available_W_per_pulse_A_B=WATT,source_available_W_per_pulse_C=1/16,
           exact_raw_identity="d_C^2-d_B^2 = (1/8) sum_over_non_B_offdiagonal_edges |Delta_S|^2 + (1/16) sum_diagonal |Delta_S|^2",
           deterministic_ball_projection="P=(C_rt+C_tr)/sqrt(2) is a contraction on channel noise; the declared C radii are twice B, so this yields only a 2r_B ball guarantee",
           independent_bare_cost=bare_cost,full_timing=full["timing"], power=full["power"])
    for key in ["total_transpose_reciprocity_rel","calibrated_transpose_reciprocity_rel","A_role_reversal_rel",
                "A_direct_vs_projection_rel","bare_VIE_loaded_total_rel","bare_VIE_loaded_calibrated_rel",
                "gaussian_C_symmetrized_signal_equals_B_rel"]:
        assert z[key] < 1e-8, (key,z[key])
    assert len(TX)==6435 and z["A_B_reciprocal_overlap"]==32
    return z,full


def summarize_heldout(powers, metrics):
    rows = []
    truth = powers["mixed_exact"]
    names = list(powers) + list(metrics)
    for h, case in enumerate(CASES):
        ix=np.arange(len(CASES))!=h
        def select(score):
            return int(np.argmax(score))
        values={"A":IA,"B":IB,"uniform_angles":IB,
                "random":int(np.random.default_rng(20261128+h).integers(len(TX))),
                "heldout_pair_ceiling":select(truth[h])}
        for name in names:
            pool=(powers if name in powers else metrics)[name][ix]
            # Train and test on independent whole material objects. Equal weight across
            # training objects by rank; raw powers are used for explicit maximin control.
            rank=np.array([rankdata(r,method="average") for r in pool])
            values[name+"_median_rank"]=select(np.median(rank,axis=0))
            if name in powers:
                values[name+"_maximin"]=select(np.min(pool,axis=0))
        r=rankdata(truth[h],method="max")
        for method, a in values.items():
            rows.append(dict(case=case,method=method,cut=a,heldout_D_over_2r=float(np.sqrt(truth[h,a])),
                heldout_fraction_of_pair_ceiling=float(np.sqrt(truth[h,a]/np.max(truth[h]))),
                heldout_worse_cut_fraction=float((r[a]-1)/(len(TX)-1))))
    return rows


def run():
    start=now()
    engine=PortEngine()
    geometry,q=bistatic_geometry(engine.op["positions"],float(engine.op["k"]))
    save(OUT/"geometry.npz",positions=engine.op["positions"],tx=TX,rx=RX,edge_i=EDGE_I,
         edge_j=EDGE_J,cut_edges=INCIDENCE.astype(bool),bistatic_q=q,geometry_score=geometry,
         B_unit=engine.b,C_unit=engine.c,empty_S=engine.empty)
    keys=["mixed_exact","born","centered_J","hybrid","centered_secant"]
    powers={key:[] for key in keys}
    metrics={k:[] for k in ("geometry","local_fisher","nuisance_fisher","conditional_fisher","signal")}
    rows=[]
    for case in CASES:
        s=now()
        d=load(PRIOR/"results"/"aperture_control"/f"{case}.npz",
               ("zero","patch_candidate","roi","patch","A_one_way_Y","A_one_way_r"))
        p,z=d["patch_candidate"],d["zero"]
        roi=d["roi"].ravel().astype(bool)
        radius=float(d["A_one_way_r"])
        edit=p-z
        mid=(p+z)/2
        f0,f1=engine.state(z),engine.state(p)
        cm=engine.state(mid,jacobian=True)
        unit=engine.c
        b=engine.b
        rf=engine.reciprocal_scale()
        mats=dict(
            mixed_exact=rf*(f1["field"].T@(edit[:,None]*f0["field"])),
            born=unit@(edit[:,None]*b),
            centered_J=cm["jacobian"]@edit,
            hybrid=.5*(unit@(edit[:,None]*cm["field"])+rf*(cm["field"].T@(edit[:,None]*b))),
            centered_secant=.5*((f0["delta"]-cm["delta"])*-2+(f1["delta"]-cm["delta"])*2))
        replay=f1["delta"]-f0["delta"]
        identity=float(np.linalg.norm(replay-mats["mixed_exact"])/np.linalg.norm(replay))
        assert identity<1e-9
        for key, matrix in mats.items():
            powers[key].append(cut_power(matrix)/(2*radius)**2)
        lf,nf=fisher_scores(cm["jacobian"],engine.op["centers"],roi,radius)
        metrics["local_fisher"].append(lf)
        metrics["nuisance_fisher"].append(nf)
        metrics["conditional_fisher"].append(exact_conditional_fisher(cm["jacobian"],engine.op["centers"],roi,radius))
        metrics["geometry"].append(geometry)
        metrics["signal"].append(cut_power(cm["delta"]))
        ep=abs(edges(replay))**2*WATT
        a,bcut=INCIDENCE[[IA,IB]].astype(bool)
        same_a_side=~a
        c_raw=float(np.linalg.norm(replay[:16,:16]/4)**2)
        b_raw=float(ep[bcut].sum())
        c_identity=float(ep[~bcut].sum()+np.sum(abs(np.diag(replay)[:16])**2)/16)
        assert c_raw >= b_raw and abs((c_raw-b_raw)-c_identity)<1e-12
        oracle_blocks=replay[RX[:,:,None],TX[:,None,:]]
        sigma=np.linalg.svd(oracle_blocks,compute_uv=False)[:,0]
        oracle=int(np.argmax(sigma))
        _,_,vh=np.linalg.svd(oracle_blocks[oracle],full_matrices=False)
        u=vh[0].conj()
        measured=(oracle_blocks[oracle]@u)
        assert abs(np.linalg.norm(measured)-sigma[oracle])<1e-12
        a_power=float(ep[a].sum())
        parts=dict(shared_A_B=float(ep[a&bcut].sum()),new_B_same_A_side=float(ep[same_a_side&bcut].sum()),
              removed_A_cross_side=float(ep[a&~bcut].sum()),new_B_cross_A_side=0.0,
              new_B_same_A_tx_arc=float(ep[~a&bcut&(EDGE_J<8)].sum()),
              new_B_same_A_rx_arc=float(ep[~a&bcut&(EDGE_I>=8)].sum()))
        labels=dict(A=np.sqrt(a_power)/(2*radius), B=np.sqrt(b_raw)/(2*radius),
                    C=np.sqrt(c_raw)/(4*radius))
        row=dict(case=case,observation_id=f"bmtm26_{case}_f400_a08_c48",
             native_current_defects=[f0["state_defect"],f1["state_defect"],cm["state_defect"]],
             mixed_identity_relative_error=identity, roi_cells=int(roi.sum()),connected_patch_cells=len(d["patch"]),
             ROI_mean_pair_difference=float(edit[roi].mean()),radius=radius,
             fixed_ratios={key:float(val) for key,val in labels.items()},energy_A_raw=a_power,
             energy_B_raw=b_raw,energy_C_raw=c_raw,energy_C_minus_B_identity=c_identity,edge_partition=parts,
             B_new_fraction=float(parts["new_B_same_A_side"]/b_raw),
             best_cut=int(np.argmax(powers["mixed_exact"][-1])),
             best_ratio=float(np.sqrt(np.max(powers["mixed_exact"][-1]))),
             coherent_oracle=dict(cut=oracle,pulses=1,available_W=1,complex_reads=8,disjoint_ratio_if_same_absolute_per_read_ball=float(sigma[oracle]/(2*radius/np.sqrt(8))),
                 ratio_at_original_64_read_radius=float(sigma[oracle]/(2*radius))),
             metrics_relative_to_mixed={
                 key:dict(complex_rel=float(np.linalg.norm(matrix-replay)/np.linalg.norm(replay)),
                     all_cut_rank_rho=float(spearmanr(powers[key][-1],powers["mixed_exact"][-1]).statistic))
                 for key,matrix in mats.items() if key!="mixed_exact"},**elapsed(s))
        rows.append(row)
        save(OUT/f"{case}.npz",exact_replay=replay,original_zero=z,original_patch=p,roi=roi,
             patch=d["patch"],original_noisy_A=d["A_one_way_Y"],endpoint_zero=f0["delta"],
             endpoint_patch=f1["delta"],coherent_input=u,coherent_cut=oracle,
             coherent_output_difference=measured,coherent_sigma_by_cut=sigma,**mats)
        if case==CASES[0]:
            receipt, full = fairness(engine,z,replay,radius,np.random.default_rng(5812026))
            dump(OUT/"port_fairness.json",receipt)
            save(OUT/"port_fairness_raw.npz",total=full["total"],delta=full["delta"],empty=full["empty"],
                 current=full["current"],loaded_antenna=full["antenna"])
        print(case,row["fixed_ratios"],row["best_cut"],row["best_ratio"],row["B_new_fraction"],flush=True)
    powers={k:np.asarray(v) for k,v in powers.items()}
    metrics={k:np.asarray(v) for k,v in metrics.items()}
    save(OUT/"cut_scores.npz",**powers,**metrics)
    heldout=summarize_heldout(powers,metrics)
    with (OUT/"heldout.csv").open("w") as f:
        writer=csv.DictWriter(f,fieldnames=list(heldout[0]))
        writer.writeheader()
        writer.writerows(heldout)
    methods=sorted({r["method"] for r in heldout})
    summary={m:{k:float(np.median([r[k] for r in heldout if r["method"]==m])) for k in
             ("heldout_D_over_2r","heldout_fraction_of_pair_ceiling","heldout_worse_cut_fraction")} for m in methods}
    dump(OUT/"summary.json",dict(cases=rows,heldout_medians=summary,cut_count=len(TX),
                       accounting=dict(engine.counters),**elapsed(start)))


if __name__=="__main__":
    with threadpool_limits(limits=1,user_api="blas"):
        run()
