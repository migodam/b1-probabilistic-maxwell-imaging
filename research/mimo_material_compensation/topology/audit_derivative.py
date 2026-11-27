"""Independent full-wave derivative and phase-conjugate power checks."""
from ports import A, PortEngine, dump, elapsed, load, now, project
from rank import OUT
import numpy as np
from threadpoolctl import threadpool_info,threadpool_limits


def run():
    start=now()
    engine=PortEngine()
    material=load(OUT/"n32_c01.npz",("original_zero",))["original_zero"]
    rng=np.random.default_rng(582401)
    h=rng.normal(size=engine.n)
    # The derivative is mathematical at a box boundary. The centered numerical
    # comparison uses an interior material to keep both points physically legal.
    m=(.95*material+.1)
    h*=1/np.max(abs(h))
    s=engine.state(m,jacobian=True)
    j=s["jacobian"]
    matrix=j@h
    fd=[]
    for t in (1e-3,1e-4,1e-5):
        mixed=(engine.state(m+t*h)["delta"]-engine.state(m-t*h)["delta"])/(2*t)
        fd.append(dict(step=t,relative_FD_defect=float(np.linalg.norm(matrix-mixed)/np.linalg.norm(matrix))))
    v=rng.normal(size=(24,24))+1j*rng.normal(size=(24,24))
    rhs=np.vdot(v,matrix).real
    adj=np.einsum("rji,rj->i",j.conj(),v).real
    lhs=float(adj@h)
    cut=load(OUT/"n32_c01.npz",("coherent_input","coherent_cut","mixed_exact"))
    from ports import TX,RX
    tx,rx=TX[int(cut["coherent_cut"])],RX[int(cut["coherent_cut"])]
    a=np.zeros(24,complex)
    a[tx]=cut["coherent_input"]
    full=engine.full_response(m)
    reflected=np.linalg.norm(full["total"][tx]@a)**2
    dissipated_other=np.linalg.norm(full["total"][np.setdiff1d(np.arange(24),tx)]@a)**2
    accepted=1-reflected
    incident=float(np.linalg.norm(a)**2)
    assert fd[-1]["relative_FD_defect"]<2e-7 and abs(lhs-rhs)<2e-8
    dump(OUT/"derivative_and_coherent_receipt.json",dict(JVP_centered_difference=fd,
        realification_adjoint_pair=dict(material_inner_product=lhs,complex_real_inner_product=float(rhs),abs_error=float(abs(lhs-rhs))),
        coherent_mechanism=dict(simultaneously_driven_ports=tx.tolist(),received_channels=rx.tolist(),
             available_unit_pulse_energy=incident,driven_reflected_energy=float(reflected),
             accepted_energy_all_driven_ports=float(accepted),other_matched_terminations=float(dissipated_other),
             net_radiated_plus_ohmic_energy=float(accepted-dissipated_other),
             constraint="one simultaneous 8-port coherent launch at total 1 W, unlike 8 sequential 1/8 W launches; phase vector is oracle and its hardware amplitude permission is not assumed"),
        threads=threadpool_info(),**elapsed(start)))


with threadpool_limits(limits=1,user_api="blas"):
    run()
