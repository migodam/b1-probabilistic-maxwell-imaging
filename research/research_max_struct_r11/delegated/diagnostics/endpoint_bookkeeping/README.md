# 公开封存后的端点投影梯度记账

状态：helper已实现；只用已解盲16对象集的一条预定保存record核对匹配逻辑。**主48尚未运行此helper。** 无physics/optimizer导入，无GT读取、forward或adjoint，也不改runner、scorer、选根、预算或算法。

## 原始gradient单位已核对

- `research/research_fbg_r3_experiments/src/physics.py:79–102`：默认beta=None，返回观测归一化数据项的gradient，再加原prior的greg；值也是data+reg。因此保存的是原Phi未缩放、含prior的材料gradient。
- `research/research_max_struct_r6/src/benchmark_solver.py:147–152`：先将上述原gradient写入paid record。
- 同文件175–179：observer返回原value/g之后，只有给minimize的返回值乘scale=1e6；保存record没有乘1e6。
- `experiments/twod_run.py:109–115`：ACTUAL_PAID_INCUMBENT和LAST_PAID_STATE会保存其已有gradient，cap NPZ仅保存材料/current/预测，不额外计算gradient。

## 严格匹配与指标

比较最终公开`R11_CAP_PREDICTION.npz`的576维m与同一实际输出目录下ACTUAL_PAID_INCUMBENT、LAST_PAID_STATE的m，使用`np.array_equal`逐项精确比较；不使用allclose/容差/另一个近端点。优先取完全匹配且含合法原gradient的ACTUAL记录，否则再查完全匹配的LAST。两个候选的逐项原因均保留。缺记录、不匹配或没有合法gradient时，全部所需指标为null并写原因，不从其他历史点补梯度，不调用物理。

合法匹配后，g不作任何缩放，输出：

- `projected_gradient_inf = ||m-clip(m-g,.02,20)||inf`；
- `raw_gradient_inf = ||g||inf`；
- `exact_lower_active/exact_upper_active`：m严格等于.02/20的cell计数；
- `near_lower_1e_minus8/near_upper_1e_minus8`：绝对边界距离<=1e-8的cell计数，含精确活动cell；
- 错号精确活动乘子：下界g<0、上界g>0的计数及合计；不将near标记替代精确活动集。

错号计数不附自定数值阈值，极小浮点符号仍可能被计入；这些只是事后端点描述，不是新的成功/收敛资格门或严格KKT证书。保存梯度存在率及null原因必须一起报告。

## 后续主288执行方式（现在未执行）

只有主288全部公开封存后，twod owner可运行：

```sh
"/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" research/research_max_struct_r11/delegated/diagnostics/endpoint_bookkeeping/endpoint_bookkeeping.py --manifest research/research_max_struct_r11/configs/twod_independent48_all_tasks.json --public-seal research/research_max_struct_r11/delegated/twod/independent48/secondary/public_seal/RECEIPT.json --output research/research_max_struct_r11/delegated/diagnostics/endpoint_bookkeeping/main48_after_public_seal
```

主入口要求原288唯一task IDs/48对象，并核验已有all288公开seal/gate元数据后才读取端点数组。它复用`experiments/twod_recovery.py::active_task`的纯metadata active-attempt覆盖，自动转向真实成功retry输出目录，保留原scientific identity/原目录及engineering_attempt_id，绝不回用已中断目录。Windows分隔符会规范化。

输出ENDPOINT_BOOKKEEPING.json/CSV和receipt。外部完整childCPU应由调用owner再计一次；其中自计时为嵌套信息。后续按task_id/case_id/method关联已授权GT评分表，只作事后描述，不能用这些梯度挑候选、改停止规则或修改主表成功门。

## 唯一逻辑核对

`validate_one_record.py`只用旧16公开seal的第一条预定记录`r11_bothblocks16_0000__cie_single__480CPU`，没有搜索有利对象。匹配ACTUAL/LAST均成功；投影梯度6.4221161e-7，rawgradient7.4178666e-6，精确下/上界306/154、near1e-8同数，错号精确活动乘子1+5=6。只用同一记录作one-ULP不匹配及缺gradient拒绝控制，均PASS，0物理调用/0GT。外部完整CPU0.111951s；见[一次核对回执](ONE_RECORD_LOGIC_CHECK.json)与[外部费用](VALIDATION_EXTERNAL_COST.json)。这不代表主48结果或新的机制实验。
