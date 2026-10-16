# L/Z/X 附件阅读审计（2026-10-16）

本分工完整阅读三个包的全部文本、数学、代码、schema与既有结果，共 **41个文件**。只读工作；没有运行包内检查脚本，没有Maxwell/Mie/JVP/VJP/adjoint/factorization，没有GPU、训练或网络检索。未启动子agent。本文件记录阅读，不将附件中的执行命令当作用户已授权的新实验。

包根：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research_inputs/lxz_20261016/`。

## 正式L报告确实存在

当前L包含 `L/Research_L/RESEARCH_L_FULL_REPORT_ZH.md`（65,730 bytes）和同名PDF，以及proof/experiment/literature附录、代码和接口。**本轮完整读了正式Markdown全文及全部附录。** PDF按重复正文跳过；这不冒称对PDF做版式或逐页验证。

X包在2026-10-15收到的是另一份四文件 `Research_L_Algebra_and_Interface(1).zip`。所以X的“缺正式L报告”是其旧输入事实，**不适用于本轮新L包**。旧X原文和manifest应保留，本轮审计更正即可，不改写历史。

## 逐文件阅读记录

下表路径均相对上面的附件包根；41项均已完整阅读。大报告按连续段读取，工具最初显示截断的L/X段已另行补读。

| 文件 | bytes | 行数 | 阅读内容与状态 |
|---|---:|---:|---|
| `L/Research_L/README.md` | 3437 | 41 | 包入口/状态说明 |
| `L/Research_L/RESEARCH_L_FULL_REPORT_ZH.md` | 65730 | 1015 | 正式L全文；五种局部对象、L1–L3、三算法、反例与证据边界 |
| `L/Research_L/appendix/EXPERIMENT_SPEC_ZH.md` | 4318 | 57 | L未执行E0–E3及未来建议预算，不构成授权 |
| `L/Research_L/appendix/LITERATURE_AND_EVIDENCE.md` | 4521 | 36 | L来源及术语/证据限制 |
| `L/Research_L/appendix/THEORY_APPENDIX_ZH.md` | 5898 | 110 | L条件命题证明与四个强制反例 |
| `L/Research_L/appendix/build_pdf.sh` | 545 | 15 | 仅文档编译脚本，本轮未执行 |
| `L/Research_L/appendix/latex_header.tex` | 451 | 16 | PDF排版配置，不影响物理算法 |
| `L/Research_L/code/algebra_checks.py` | 8183 | 166 | 合成Schur/precision/finite-feedback/coarea检查实现；本轮只读 |
| `L/Research_L/code/algebra_results.json` | 3302 | 119 | L既有合成结果；不是Maxwell数据 |
| `L/Research_L/interface/LocalFiberPackage.schema.json` | 11427 | 537 | research-l-0.1单独sidecar；未附真实instance |
| `L/Research_L/interface/PRO_Z_HANDOFF.md` | 3882 | 40 | L→Z/X责任与测度交接语义 |
| `Z/Research_Z_Physics_Context_Latent_Generative_Maxwell_Imaging/CODEX_NEXT_ROUND_ZH.md` | 1484 | 9 | 附件内下一阶段建议/约束；不是自动扩大当前用户授权 |
| `Z/Research_Z_Physics_Context_Latent_Generative_Maxwell_Imaging/README.md` | 3235 | 31 | 包入口/状态说明 |
| `Z/Research_Z_Physics_Context_Latent_Generative_Maxwell_Imaging/RESEARCH_Z_REPORT_ZH.md` | 25559 | 321 | 正式Z全文；prior/context/proposal分离、PZ新支持、global/local路线 |
| `Z/Research_Z_Physics_Context_Latent_Generative_Maxwell_Imaging/algorithms/ALGORITHM_SPEC_ZH.md` | 5357 | 73 | Z-PCS伪代码/spec；非runnable Maxwell实现 |
| `Z/Research_Z_Physics_Context_Latent_Generative_Maxwell_Imaging/experiments/Z0_Z3_FALSIFIABLE_PROTOCOL_ZH.md` | 7426 | 101 | Z全部DESIGN_ONLY；Z0材料表达与Z1–Z3前置条件 |
| `Z/Research_Z_Physics_Context_Latent_Generative_Maxwell_Imaging/interfaces/PRO_L_CANDIDATE_INTERFACE.schema.json` | 4580 | 96 | ZL-0.1候选schema；非L schema的同名替代 |
| `Z/Research_Z_Physics_Context_Latent_Generative_Maxwell_Imaging/literature/LITERATURE_AND_NOVELTY_AUDIT_ZH.md` | 5757 | 50 | Z近邻与不能主张新颖的通用组合 |
| `Z/Research_Z_Physics_Context_Latent_Generative_Maxwell_Imaging/provenance/CLAIM_LEDGER.md` | 3479 | 30 | Z历史/数学/假说分列 |
| `Z/Research_Z_Physics_Context_Latent_Generative_Maxwell_Imaging/provenance/SOURCE_MANIFEST.json` | 4743 | 100 | Z可变main来源指针；不是实际固定commit的完整外部证据 |
| `Z/Research_Z_Physics_Context_Latent_Generative_Maxwell_Imaging/theory/THEORY_NOTES_ZH.md` | 6518 | 95 | Z标准事实、条件推导与明确未测假设 |
| `X/Research_X/CODEX_NEXT_ROUND_ZH.md` | 6068 | 58 | 附件内下一阶段建议/约束；不是自动扩大当前用户授权 |
| `X/Research_X/README_ZH.md` | 2287 | 33 | X包入口/状态说明 |
| `X/Research_X/RESEARCH_X_REPORT_ZH.md` | 51427 | 664 | 正式X全文；冻结local/global核与固定材料current-refinement主张 |
| `X/Research_X/checks/L_ALGEBRA_RERUN.json` | 3302 | 119 | 既有L合成重放结果，不是新独立Maxwell结果 |
| `X/Research_X/checks/L_ALGEBRA_RERUN.stdout.txt` | 3303 | 119 | 同次L重放stdout，不能当独立重复 |
| `X/Research_X/checks/X_SYNTHETIC_CHECKS.json` | 1828 | 53 | 既有X合成结果；oracle toy区间使用toy真值不能作省算证据 |
| `X/Research_X/checks/X_SYNTHETIC_CHECKS.stdout.txt` | 1829 | 53 | 同次X结果stdout，不能当独立重复 |
| `X/Research_X/code/ADAPTER_CONTRACT_ZH.md` | 3120 | 39 | uniform→Gaussian prior transport及密度Jacobian合同 |
| `X/Research_X/code/requirements.txt` | 6 | 1 | NumPy通用核依赖；不代表已具备physics后端 |
| `X/Research_X/code/run_synthetic_checks.py` | 7172 | 161 | X合成/五状态/1200步toy测试实现；本轮未运行 |
| `X/Research_X/code/sampler_core.py` | 7805 | 187 | 实际冻结Gaussian-prior MH/pCN/selector/严格区间拒绝通用核 |
| `X/Research_X/experiments/DECISIVE_PROTOCOL_ZH.md` | 11485 | 132 | X E0–E4；当前新physics/adjoint/factorization/GPU均0 |
| `X/Research_X/literature/NEAREST_PRECEDENTS_ZH.md` | 18906 | 196 | 28项定向文献审计；区分摘要与正文证据深度 |
| `X/Research_X/provenance/CLAIM_LEDGER.json` | 2535 | 84 | X逐项claim等级与未证声明 |
| `X/Research_X/provenance/DELIVERY_MANIFEST.json` | 3131 | 100 | X历史交付清单；既有hash字段仅作为源文本读取，未核对/重算 |
| `X/Research_X/provenance/INPUT_MANIFEST.json` | 3621 | 89 | X所接收旧L supplement的历史清单；本轮L版本已有不同 |
| `X/Research_X/provenance/RESOURCE_LEDGER.json` | 373 | 12 | X包自报新physics0、仅synthetic；不含当前执行新费用 |
| `X/Research_X/provenance/SOURCES.json` | 8606 | 219 | X来源和实际查阅深度声明；本轮未联网重新验证 |
| `X/Research_X/provenance/SOURCE_INDEX.md` | 6242 | 213 | X原文入口；本轮只读索引，不声称读了外部全文 |
| `X/Research_X/theory/FOUR_PROPOSITIONS_ZH.md` | 15887 | 331 | X四命题、假设、证明、反例与standard/conditional区分 |

## 有意跳过的两个文件

- `L/Research_L/RESEARCH_L_FULL_REPORT_ZH.pdf`：与MD重复正文，本轮使用完整MD；未作PDF视觉QA。
- `L/Research_L/FILE_MANIFEST_SHA256.json`：按用户偏好与父任务跳过。没有计算或核对任何SHA-256。其他provenance里既有hash字段仅作为原始文档文字读取，不作为通过条件，也未重新计算。

## 可复用实现与证据强度

L仅有合成代数核验；Z没有sampler/decoder训练实现；X有通用冻结MH/pCN代码。三个包都没有实际新的Maxwell候选数组、真实LocalFiberPackage instance、可认证forward误差界或新posterior reference。既有JSON/stdout副本不是独立实验。

外部论文证据深度以原包自报为准。本次没有访问任何外部链接，不能把原包的摘要检索升级为本轮完整全文核验；尤其X的R10/DDMC和R28仅摘要，P6在SOURCE_INDEX与SOURCES中对查阅深度描述也不同，宜用较保守的“关键段/入口”而非完整全文。

父线程保留核心数学、target选择、实验门槛和最终科学裁决。
