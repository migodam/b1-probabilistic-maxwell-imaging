# R9经典二维实现核对（不使用GT）

核对对象：`src/classic_extension.py`、`experiments/classic_extension.py`、冻结的 `configs/classic_extension_main.json`。这里只核对实现和合同，不判断新理论成立或创新性。真正运行状态以queue和逐项回执为准。

| 项 | 核对结果与证据 |
|---|---|
| 原观测复用 | freeze只copy2 R8的PUBLIC_DATA和DATA_RECEIPT；solve加载R9相对路径副本，验证40×20、forward_side150和400MHz。旧输入inverse_side32仅为旧build元数据；新acquisition显式side64，不会被该旧值降到32² |
| 材料自由度 | `m0=np.zeros(model.n)`，64²全4096个χ，L-BFGS-B逐像素相同[0,20]；没有support mask、已知分区或几何参数化 |
| GlobalDual／A／A+的current | 三者传入同一固定first β、`sourcewise_cg=False`。R8 `run_cie_fft`把current_transform固定为scalar；A/A+仅material_transform变onehop。后续m不同自然使current迭代值不同，不能要求它们整条current轨迹相同 |
| Single／Dual | sourcewise PR与每源line step；single βcurrent随βmaterial下降，Dual βcurrent固定首轮β。global用整体PR/line step及Frobenius系数停止，区别保留 |
| source权重 | 未覆盖R8默认data 1/||Y_s||²、state 1/||E_s||²；固定真实观测，不按材料GT调权。不进行源混合或重新定义weights |
| γ传输 | 复用R8按frequency tuple传递旧系数、为新增频率补零；未对APCS重新做SVD白化，保留原native Fourier坐标及DPCS投影 |
| matched6 prior | 相同smooth=TV=1e−4，作用于原χ；A+用相同regularizer对象作MM接受。未把材料prior移到modified contrast |
| independent NIE／LSIE | 独立sourcewise PR闭式real modified-contrast LS；LSIE残差直接J−χQ，不是β重命名。prior=False；不能与matched6声称仅方程不同，更不是作者软件复现 |
| 全物理corrector | 所有臂同一300 accepted ceiling、1200函数评估上限、相同容差/线搜索；80检查点取同一次连续L-BFGS-B，无重启。最终保留已付原物理Phi最低incumbent，同时保存optimizer终点用于核查 |
| exact LU fallback | 保留R8 GMRES≤1000、rtol1e−9和实际残差≤2e−9检查；缓存的是不变G，仍构造同一I−diag(χ)G并库LU。source／adjoint共用同m LU；失败迭代与重复RHS均计入。无失败state的伪有限目标／零梯度 |
| 成本 | child由wait4获得全进程user+sys CPU与peak RSS；内部计时含入口后的导入/组装，供Phi曲线。账本不把aux CPU、物理组件CPU再加到process CPU；Green／coefficient FFT和RHS计数分列 |
| GT顺序 | solve/src中没有GT数组加载或构造。score先检查全部16项完成封存／明确失败，再写已付public-Phi CPU前缀选择，才构造公开Austria真值并写sealed_truth。不会读R9新确认GT |
| CPU前缀 | 300/600/900在GT前冻结。只用已有stage/aux/80/final及真实规定零初值；保守上移含启动/输出尾部，输出最新保存点与gap。未保存材料的实际incumbent未知，不能称精确等CPU截断实验 |

明确限制：paper stopping可能被500-stage cap截断，physical corrector可能300 CAP；这两种都必须保留停因，不称收敛。64²仍有约9.96%的已知400MHz clean跨网格误差；固定10%噪声不提供离散误差或材料恢复保证。所有对象属于公开论文开发案例。noise-rule rank补充为NOT_RUN，未因结果追加计算。

2019原文Eq.(16)/(17)根号抽取错误已另记在 `AUDIT.md`、`CIE2019_FORMULA_EVIDENCE.xml` 与R9数据包corrections文件；当前二维代码从始至终使用正确的二维δ²门限，未受此元数据错误影响。
