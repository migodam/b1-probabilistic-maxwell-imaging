Findings:

1. Package integrity passes: all 25 manifest-listed files match their recorded byte counts and SHA-256 hashes; all JSON parses and there are no unlisted or unsafe paths.

2. Historical commit `55724cb1a7bafbd9afa3b95ac817522f8dcaf2d0` exists locally. R2-R4 trees have no content difference from that commit to HEAD `18620d199a4882c1f189dcd276c6572f7414e353`.

3. `GENERAL6_00` observation: `research/research_z_alg_r2/experiments/general6_public/GENERAL6_00.npz`, containing only `y` with shape `(1,3)`, complex128. Geometry is in `experiments/GENERAL6_CONTRACT.json`.

4. Root provenance is exact:
   - Q20: `results/general6_GENERAL6_00_general_quotient20.json` — 20 endpoints, 2 qualified roots.
   - 64-start baseline: `results/general6_GENERAL6_00_ReducedFullJ64.json` — 1 qualified root.
   - Root separation: `1.1587903026947952` material RMS.
   - Status remains `NOT_RIGOROUSLY_CERTIFIED`.

5. Source API gap: existing code supports point-source positions and fixed Hankel fields, or precomputed `G/E0/S`. No audited API provides transmitter apertures, complex excitation control, coherent source mixtures, power/energy normalization, or a source Gram matrix. This blocks treating an arbitrary support-space `b` as physically realizable.

6. Count provenance is verified: R2 `2/1/1` versus `1/1/1`, 192 starts and 24 cap hits; R3 `864` Q6 calls, `538` new legal roots, `462` distant roots, zero direct global-best updates; R4 basin `29+18` runs, diversity `360=288+72`.

7. All 24 bibliography labels L1-L24 are present locally. No DOI, arXiv record, or URL was fetched, so external verification remains untested. L17 has no stable identifier in the package; L2 and L18 were only abstract/metadata-checked.

8. The sparse checkout omits `research/research_z_alg_asymmetric_validation`, although it exists in the reference tree/index. Separately, untracked `research/research_z_alg_r5_measurement/` files appeared during the audit and were not inspected or modified.

Artifacts:

- [audit.md](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/z_alg_r5_input_audit/audit.md)
- [evidence.csv](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/z_alg_r5_input_audit/evidence.csv)
- [bibliography.json](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/z_alg_r5_input_audit/bibliography.json)
- [SUMMARY.md](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/z_alg_r5_input_audit/SUMMARY.md)

No experiment was run and no historical file was changed.