# GN–NN R3: Voxel versus receiver-null current material completion

Completed full nonlinear scalar-TM VIE study: 1024 training mothers,64 validation mothers,64 independent test mothers /192 paired observations,17 confirmation models. Primary full-voxel/full-current/full-null comparison,three training seeds,two seed7 cascades,complete physical correction,and full-operator filter ablations.

Measured outcome: hard receiver-null projection did not improve noisy material recovery over full-current supervision. V-physics has a strong seed7 result but weaker seeds17/27. J-full was more stable across the three seeds. Two-pass voxel completion is promising at one seed. No high-contrastA12,probability/root-enumeration or universal superiority claim.

- [Chinese report](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_gn_nn_completion_r3/)
- [Raw theory handoff](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_gn_nn_completion_r3/CODEX_TO_PRO_HANDOFF.md)
- [All test images](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_gn_nn_completion_r3/figures/all_test_cases/GALLERY.html)

Download and extract all six ZIP parts into one directory,preserving research/ paths. Inputs/anchors,trained weights,predictions/scores,costs/logs,source/reports and direct frozen physics dependencies are included. RELEASE_MANIFEST.json lists sizes and scope. Reproduction does not require downloading new MNIST or replacing the full-VIE physics.

Known limitations retained: digit–amplitude correlation,training seed sensitivity,continuous-material Otsu calibration,GN budget/numerical stops,and one failed attempt with unknown LU/RHS. Truth is now publicly disclosed; historical re-scoring is not a new blind test. No hash checks were used.
