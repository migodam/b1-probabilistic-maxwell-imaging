# GN–NN R3

COMPLETED:1024training mothers,64validation mothers and64independentTEST mothers;17models, full-VIE correction and three training seeds. No Weak64/128 or hard weak-spectrum truncation.

[中文实验报告](RESEARCH_REPORT_ZH.md) · [理论交接](CODEX_TO_PRO_HANDOFF.md) · [全部测试图](figures/all_test_cases/GALLERY.html) · [复现](REPRODUCE.md) · [下载数据与权重](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/gn-nn-completion-r3-20261107)

Main finding: receiver-null current learning did not improve noisy recovery over full-current material supervision across three seeds. Voxel physics input has a strong seed7 result but training instability in seeds17/27. Two-pass voxel completion is promising at one seed. These are synthetic low-contrast full-VIE results, not a high-contrast/global-root claim.

Frozen protocol: [audit/PROTOCOL.json](audit/PROTOCOL.json). [Final scientific review](audit/FINAL_SCIENTIFIC_REVIEW.md).
