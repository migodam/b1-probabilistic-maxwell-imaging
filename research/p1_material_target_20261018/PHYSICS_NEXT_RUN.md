> 历史冻结建议：用户随后明确授权，本资格阶段现已完成并封账。实际结果为26states/976RHS、1/2合格；见[真实物理报告](PHYSICS_QUALIFICATION_REPORT_ZH.md)。以下原预算建议保留，不代表当前仍未授权，也不授予后续比较额度。

# 第一轮P1真实物理：代码已具备，预算尚未授权

本阶段只生成材料和冻结target，Maxwell新增硬上限0。EXECUTION_AUTHORIZATION保持0。PHYSICS_AUTHORIZATION_TEMPLATE内的32/1024仅为建议，explicit_user_physics_authorization=false，不能直接运行。

## 固定任务

- 原scalar TM24、8sources、32complex receivers、host .1+.04i、gain1与固定sigma，不修改噪声以制造歧义。
- 两个独立development parent：生成前确定的dev000、dev001；在线观测命名obs_A/B，不提供隐含latent/seed。
- 先各付一次truth full-forward生成观测，再从共同prior均值开始独立GN/LM，最多12个新状态/parent，所有line search与拒绝计费。
- 资格为prior合法、完整状态backward residual≤1e−10及预先固定noise-ball；报告KKT/gradient/latent logprior。此处不要求已认证MAP，也不称posterior采样。
- 任一未合格就保留失败、不用GT替代，不增加预算找漂亮对象。

## 预算请求

| 项目 | states | complex RHS |
|---|---:|---:|
| 两个真实材料观测 | 2 | 16 |
| 两个独立反演，各至多12states | 24 | 192 |
| 至多24个共享receiver-adjoint J | 0 | 768 |
| 失败/数值异常预留 | 6 | 48 |
| 总上限 | **32** | **1024** |

受控physics wall600s、NN/GPU/Mie/新3D为0。state与RHS先到哪个上限就停，失败不能刷新额度。从零生成P1prior无需这些物理额度，已经完成。

获得明确新授权后应创建绑定本target/本stage的新authorization文件，再运行`PYTHONPATH=src python tools/run_p1_physics.py --authorization <新授权路径>`。同一物理ledger禁止重复启动；恢复要保留已有费用。当前0配置执行会在读取offline truth及导入backend之前终止。

## 与A17-X的关系

资格阶段不并行跑QD/CEM/大型NN，也不消耗此前A17-X建议额度。只有合法基础材料和观测齐全后，另开24states/1024RHS/600s的共同机制比较；该预算同样未授权。exact T-matrix、普通Full-J、goal-ROM不能省略。若P1物理本来可辨识或普通方法更好，接受负结果。
