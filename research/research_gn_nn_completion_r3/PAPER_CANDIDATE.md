# Paper candidate — HOLD for a standalone methods claim

**Working title:** *Receiver-Null Current Learning versus Voxel Completion in Nonlinear Electromagnetic Imaging: Shared-Material Readout, Noise and Physical Correction*

The scientific problem is whether a large receiver-null component of the GN current error identifies a better learning representation for material recovery. The experiment shows that current-error energy, material readout accuracy and full-VIE feasibility need not agree.

## At most three defensible candidate contributions

1. A controlled comparison of full voxel, full current and complete receiver-null current outputs with the same physical input data, backbone, paired mother objects, frozen training seeds and counted full-Maxwell correction. This is comparative evidence, not an invention of current learning or nullspace networks.
2. A conditional real shared-material readout Gram and its relation to a weighted annihilator; the observed575rank and negligible GT-error energy in the excluded first-order direction help rule out a simplistic dimensional-capacity explanation. The exact containment assumptions are not proved for arbitrary Maxwell discretizations.
3. An experimentally documented mismatch between projected-current supervision and material recovery, including the clean/noisy natural-oracle contrast and the benefit/cost of retaining learned-prior selection under a physical data gate. This motivates, but does not implement, a noise-aware soft receiver release.

## Direct prior-art risk

[Si et al.,2024](https://www.mdpi.com/1424-8220/24/18/5997) already uses two-stage contrast-source learning for electromagnetic inverse scattering. [Schwab, Antholzer and Haltmeier,2018](https://arxiv.org/abs/1806.06137) provides deep nullspace learning. CSI/TSOM, variable projection, shared-material LS readout, physics correction and residualCNN architectures are established tools. The present architecture alone is insufficient novelty.

## What cannot be claimed

No general hard-null advantage, no stable V-physics superiority acrossseeds, no posterior sampling, no complete root set, no global GT convergence, no high-contrastA12recovery, no proof of universal575rank, no noise-optimal classical benchmark domination, no end-to-end GPU acceleration.

## The decisive missing experiment

Freeze independently crossed shape/font/amplitude/background data and noise-aware conventional baselines, then compare hard-null, bounded visible-current release and full current using a material-relevant loss/metric derived from the readout Gram. Keep three training seeds and full physical accounting. A robust material gain from that controlled release could support a methods paper; the current study alone is a reproducible mechanism/benchmark result.
