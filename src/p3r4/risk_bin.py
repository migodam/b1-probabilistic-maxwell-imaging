"""Revision1: exact h-mass bin moment labels for the unchanged P1 decoder.

Bin-averaging analytically Rao--Blackwellizes the randomized threshold label.
This is the classical layer-cake identity, not a new theorem or posterior law.
"""
import numpy as np
import torch
from scipy.special import ndtr
from scipy.optimize import brentq
from p1_material.prior import saturation
from .risk import LocalRisk

MOMENT_RULE='equal_h_mass_bin_integrated_joint_delta_event_v1'


def quantile_bins(k):
    """Original h maps[-1,1] onto[0,1]; interior roots preserve that decoder."""
    if k not in (8,16):raise ValueError('revision1 freezes K8/16')
    edges=np.array([-1.]+[brentq(lambda s:saturation(s)-j/k,-1,1,xtol=1e-14) for j in range(1,k)]+[1.])
    centers=np.array([brentq(lambda s:saturation(s)-(j+.5)/k,-1,1,xtol=1e-14) for j in range(k)])
    return edges,centers


def joint_bin_labels_from_s(s,delta,k):
    """Analytic E_T[delta*1{s>T}|T in bin_j],T density h'.

    s[N,ROI] standardized field; delta[N,2] real/imag joint attributes.
    The same parent's delta and s are retained together in every label.
    """
    s=np.asarray(s);delta=np.asarray(delta)
    if s.ndim!=2 or delta.shape!=(len(s),2):raise ValueError('s[N,ROI],delta[N,2] required')
    h=saturation(s)  # Source decoder, including exact plateaus.
    labels=delta[:,:,None,None]*np.clip(k*h[:,None,None,:]-np.arange(k)[None,None,:,None],0,1)
    return labels


def bin_phase_labels(x,geometry,thresholds,s=None):
    """Training-only labels; thresholds selects K, not point-event locations."""
    k=len(thresholds);cfg=geometry['config'];roi=geometry['roi']
    if k not in (8,16):raise ValueError('revision1 K8/16 only')
    if s is None:s=(cfg['field_mean']+np.asarray(x[:,:64])@np.asarray(geometry['A'])[roi].T)/cfg['phase_transition_halfwidth']
    if np.asarray(s).shape!=(len(x),len(roi)):raise ValueError('standardized ROI s shape')
    u=ndtr(x[:,64:]);r=cfg['inclusion_increment_real_interval'];i=cfg['inclusion_increment_imag_interval']
    delta=np.stack((r[0]+np.diff(r)[0]*u[:,0],i[0]+np.diff(i)[0]*u[:,1]),1)
    return joint_bin_labels_from_s(s,delta,k).astype('float32')


class BinRisk(LocalRisk):
    """Same6114/6386 parameters and monotone head; only moment rule changes.

    State keys match LocalRisk exactly. RiskActionPredictor loads saved weights
    and bin-quantile centers through existing buffers and replays without edits.
    """
    def __init__(self,kind,geometry):
        if kind not in ('phase8','phase16'):raise ValueError('phase8/phase16 only')
        super().__init__(kind,geometry)
        edges,centers=quantile_bins(self.k)
        self.weights.fill_(1/self.k)
        self.thresholds.copy_(torch.tensor(centers,dtype=self.thresholds.dtype))
        self.bin_edges=edges  # Metadata, deliberately not an extra state buffer.
