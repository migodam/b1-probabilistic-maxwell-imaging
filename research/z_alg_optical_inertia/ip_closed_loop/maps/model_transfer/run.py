"""Paired discrete-grid transfer on saved 96² materials."""
import json
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

import numpy as np
from scipy.special import hankel1, jv
from threadpoolctl import threadpool_limits

from common import FREEZE, HERE, Geometry, clock, elapsed, geom_path, save_json


OUT=Path(__file__).resolve().parent
RECEIVER_FACTOR=np.sqrt(20/32)


def unit_check(fine,coarse):
    cfg=FREEZE["model"]
    apt=FREEZE["aperture"]
    k=2*np.pi*cfg["frequency_hz"]/cfg["host_speed_m_per_s"]
    norm=np.linalg.norm
    rows={}
    receiver_degrees=np.deg2rad(apt["rx_angles_degrees"])
    rx=apt["sensor_radius_m"]*np.c_[np.cos(receiver_degrees),np.sin(receiver_degrees)]
    coefficients={}
    for g in (coarse,fine):
        dx=cfg["side_m"]/g.side
        a=dx/np.sqrt(np.pi)
        xy=(g.xy+.5)*dx-cfg["side_m"]/2
        base=.5j*np.pi*k*a*jv(1,k*a)
        coeff=base*(RECEIVER_FACTOR if g.side==96 else 1.)
        coefficients[g.side]=coeff
        rho=np.linalg.norm(rx[:,None,:]-xy,axis=2)
        sensor=coeff*hankel1(0,k*rho)
        z=g.B.reshape(g.side,g.side,-1)
        xx=np.angle((z[:,1:]*z[:,:-1].conj()).mean(axis=(0,1)))/dx
        yy=np.angle((z[1:]*z[:-1].conj()).mean(axis=(0,1)))/dx
        rows[str(g.side)]=dict(incident_angles_degrees=np.rad2deg(np.arctan2(yy,xx)).tolist(),
                              incident_wavenumber_rad_per_m=np.hypot(xx,yy).tolist(),
                              transmitter_modulus=[float(abs(z).min()),float(abs(z).max())],
                              sensor_row_2norm=norm(g.S,axis=1).tolist(),
                              equal_disk_coefficient_im=float(coeff.imag),
                              sensor_array_analytic_relative_error=float(norm(g.S-sensor)/norm(sensor)))
    pos=fine.B.reshape(96,96,7)[1::3,1::3].reshape(coarse.B.shape)
    spotted=fine.S.reshape(9,96,96)[:,1::3,1::3].reshape(coarse.S.shape)
    amp=coefficients[32]/coefficients[96]
    return dict(receiver_radius_m=apt["sensor_radius_m"],tx_indices=apt["tx_indices"],
                rx_indices=apt["rx_indices"],native_receiver_quadrature_factor=float(RECEIVER_FACTOR),
                coarse_to_native_readout_factor=float(RECEIVER_FACTOR),
                coincident_centre_B_relative_error=float(norm(pos-coarse.B)/norm(coarse.B)),
                coincident_centre_S_ratio=float(amp.real),
                coincident_centre_S_scaled_relative_error=float(norm(amp*spotted-coarse.S)/norm(coarse.S)),
                port_rows=rows)


def main():
    started=clock()
    fine,coarse=Geometry(geom_path(96)),Geometry(geom_path(32))
    check=unit_check(fine,coarse)
    assert check["coincident_centre_B_relative_error"]<2e-13
    assert check["coincident_centre_S_scaled_relative_error"]<2e-13
    for x in check["port_rows"].values():
        assert x["sensor_array_analytic_relative_error"]<2e-13
    norm=np.linalg.norm
    rows=[]
    for name in ("M96-R01","M96-P01"):
        oid=f"ipcl-fov-f1__{name}__weak__7tx9rx__n03"
        start=clock()
        ref=np.load(HERE/"offline"/f"{oid}.npz",allow_pickle=False)
        obs=np.load(HERE/"public"/f"{oid}.npz",allow_pickle=False)
        x96=ref["truth"]
        x32=x96.reshape(32,3,32,3).mean(axis=(1,3)).ravel()
        ref_forward,j,info=fine.forward(x96)
        assert info==0
        assert norm(ref_forward-ref["clean"])<1e-10
        Y32,j32,info=coarse.forward(x32)
        assert info==0
        comparable=RECEIVER_FACTOR*Y32
        Y96=ref["clean"]
        noise=float(obs["r"])
        mismatch=float(norm(comparable-Y96))
        residual=float(norm(comparable-obs["Y"]))
        raw_difference=float(norm(Y32-Y96))
        r32=noise+mismatch
        trace96=j-x96[:,None]*(fine.B+fine.apply(j))
        trace32=j32-x32[:,None]*(coarse.B+coarse.apply(j32))
        np.savez_compressed(OUT/f"{name}.npz",mat32=x32.reshape(32,32),fine_clean=Y96,
                            coarse_clean_native_units=comparable,coarse_clean_unweighted=Y32,
                            original_observation=obs["Y"],fine_contrast=x96.reshape(96,96),fine_noise_radius=noise,
                            coarse_forward_residual_norm=norm(trace32))
        row=dict(object_id=name,observation_id=oid,material_interpolation="3×3 box area-average",
                 measurement_norm=float(norm(Y96)),original_noise_radius=noise,
                 coarse_minus_native_frobenius=mismatch,
                 relative_model_error_vs_clean=mismatch/float(norm(Y96)),
                 model_difference_over_noise_radius=mismatch/noise,
                 unharmonized_difference_over_noise_radius=raw_difference/noise,
                 direct_coarse_data_residual=residual,direct_coarse_data_residual_over_noise=residual/noise,
                 sufficient_paired_radius=r32,sufficient_paired_radius_over_native=r32/noise,
                 coarse_state_residual_fro=float(norm(trace32)),fine_state_residual_fro=float(norm(trace96)),
                 **elapsed(start))
        rows.append(row)
        print(name,row,flush=True)
    save_json(OUT/"transfer.json",dict(scope="Paired actual materials; no uniform bound over either material class",
                                      unit_check=check,paired_results=rows,timing=elapsed(started)))


if __name__=="__main__":
    with threadpool_limits(limits=1,user_api="blas"):
        main()
