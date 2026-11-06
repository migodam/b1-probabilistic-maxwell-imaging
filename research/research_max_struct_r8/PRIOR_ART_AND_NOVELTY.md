# R8 最近先例与贡献边界（2026-11-06）

当前证据不支持把“NN 给初值、传统全波反演负责最终收敛”或“先学习辅助电流、再进行 TSOM/CSI”作为 R8 的新概念。R8 值得检验的窄贡献是：给定独立对象/物理条件划分与完整训练成本后，轻量材料初始化是否与特定 A+ 材料读出产生可重复、可解释的计算或恢复收益。是否成立由本轮封存后实验决定；此处只是先例边界，最终科学判断归主线程。

## 检索与证据范围

实际按 `scholarqa-research` 运行两组互补查询，问题分别覆盖电磁 NN 初始化与 IMSA-CIE。查询的四个子请求都返回 HTTP429，空结果不能解释为没有先例。随后依用户授权读取作者/出版社原文，保留 `delegated/paper_baselines/SCHOLARQA_COLLECT.json`。最终独立运行 Semantic Scholar 批量身份验证，五篇全部 resolved，见 `SCHOLARQA_VERIFY.json`；这次身份验证成功不应被误写成原检索成功。搜索不是系统综述，也不保证全局新颖性。

| 有限判断 | 来源、证据层级、具体位置 | 对 R8 的影响 |
|---|---|---|
| NIE/CIE 变换、FFT-TSOM current regularization 和多轮 beta/Fourier continuation 已有 | Zhong 等，2016，Tier A，作者完整论文 Eq3、§III、§IV；[DOI](https://doi.org/10.1109/TAP.2016.2535492) | 不能把标量分母、beta continuation 或 Fourier 电流分解当成新家族 |
| current/material 使用不同 beta 已有，而且过大的 current beta 也有失败实例 | Zhong & Xu，2019，Tier A，§2.4、§3.2；[DOI](https://doi.org/10.3390/jimaging5020027) | Dual-CIE 是必要基线；参数放松不能全部计入 A+ 创新 |
| IMSA 与 CIE-I 的多分辨率集成已有数值及实验测试 | Zhong、Salucci、Xu、Polo、Massa，2020，Tier B，作者机构摘要；[机构记录](https://iris.unitn.it/handle/11572/260153)、[DOI](https://doi.org/10.1109/TMTT.2019.2956939) | 原创性不能建立在“把 CIE 做成多分辨率”上。本轮没有取得其受限全文，不虚构其具体收缩域/停止细节，也没有冒名执行 IMSA |
| 学习 radiation noise-space contrast source，再以 SOM/TSOM+TV refinement 已有 | Sanghvi、Kalepu、Khankhoje，IEEE TCI 6:46–56，2020（early-online DOI2019），Tier A，作者 PDF pp1–4，Fig1、§III.B–D、§IV；[作者论文](https://ee.iitm.ac.in/~uday/pub/2019-tci-dlsom.pdf)、[DOI](https://doi.org/10.1109/TCI.2019.2915580) | 这是比泛泛的“learning warm start”更近的直接先例；训练 GT 电流、合法 DPCS/APCS 初始化和 source geometry 变化均应逐项比较 |
| NN 估计初值再进行 Gauss–Newton 的 warm start 在逆声散射已有 | Zhou、Han、Rachh、Borges，2023，Tier A，作者全文 §1、§4–5；[arXiv](https://arxiv.org/html/2212.08736v3)、[DOI](https://doi.org/10.1016/j.jcp.2023.112341) | 能否改变初始化依赖、训练样本随难度增长、噪声/孔径迁移都是已有问题；其 sound-soft star-shaped 几何参数并非本轮全像素介质 Maxwell |
| U-Net 产生介电常数初值，再用原 DBIM refinement 的电磁方案已有 | Magdum、Beerappa、Erramshetty，2023 online /2024 issue，Tier B，出版社摘要；[DOI](https://doi.org/10.1515/freq-2023-0074) | 即使仅预测材料而不预测辅助电流，也不足以主张新范式。出版社网页与 PDF 多次405，已读取可检索出版社摘要，未假称取得完整网络训练细节 |
| 神经化迭代电磁逆散射也已有 | Shan 等，2023，Tier B，作者摘要，[arXiv](https://arxiv.org/abs/2112.09831)、[DOI](https://doi.org/10.1109/TAP.2022.3217333) | NeuralBIM 学习迭代更新；R8 的小模型初始化＋解析求解器需以实验区别于这类更重 unrolling，不直接声称统一优越 |
| CIE框架自身也已有小CNN＋least-squares模块的神经unrolling | Xu 等，2023，Tier B，原论文摘要；[DOI](https://doi.org/10.1109/TAP.2023.3305834)、[作者报告及引用](https://author.piers.org/ac_api/preview.php?id=231113033645&t=ab) | CIE-SOM-NET在limited aperture下学习current/modified contrast更新；泛称“NN与CIE结合”不能作为新贡献。R8是材料warm-start后保留解析current/readout，属于不同具体接口，收益仍要实测。未取得原全文/软件，不虚构其训练规模或成熟复现资格 |
| Fourier展开contrast source＋CIE＋低频到高频多轮优化已有 | Xu、Zhang、Wei，2020，Tier B，原论文摘要与作者机构列表；[DOI](https://doi.org/10.1109/TMTT.2020.2978165)、[作者机构列表](https://faculty.hdu.edu.cn/dzxxxy/xkw_en/main.htm) | 不能仅凭native FFT、Fourier mask continuation或小current subspace声称新的科学表示；本轮nativeFFT是继承工程实现 |
| 大beta CIE估计current/total field，再作regularized线性contrast恢复已有 | Bevacqua、Isernia、Crocco，2026，Tier A，本轮新取得作者完整论文§III StepsI–II Eq4–6；[DOI](https://doi.org/10.1109/LAWP.2026.3660123)、[作者全文](https://www.researchgate.net/publication/400394332_A_Low-Nonlinearity_Approach_to_Inverse_Scattering_via_an_Unconventional_Exploitation_of_the_Contraction_Integral_Equation) | LUCIE直接限制泛称“先用CIE估电流再改材料读出”的原创性；其数据到contrast线性系统与本轮固定current的Möbius条件材料块不是同一更新。全文有可选support限制；本轮基准所有像素自由，不能偷换成使用真support的比较 |

## 本轮可写、不可写的内容

可以准确写：本轮独立测试真实论文的物理尺寸、20 平面波、40 接收点和对应频率；执行独立 NIE/LSIE 方程重实现；将 A/A+ 与 current-sourcewise/global 进行配对；报告 full-wave LU/RHS 或 GMRES/FFT 工作量及完整训练成本。原论文无源码不是不能复现的理由，因此本轮已编写并实际执行独立方程实现。

不能把该新实现直接冠名为“已认证的成熟原作者 NIE 软件”。本轮明确使用实非负全像素材料、有限迭代上限、指定圆盘近似单元与部分缩小反演网格。自有 A+/FFT 接口的程序可运行，只说明实现状态；并不自动证明连续 Maxwell 精度、跨 basin 改善或论文创新。

不需要用 NN 使初始原物理目标优于均匀初值才允许后续迭代；必须实际跑共同预算。反之，得到较低残差也不能替代材料恢复评价。不同对象的成功率与配对费用、旧暴露对象/新确认对象身份，以及训练/验证/测试家族划分，应从真实实验记录读取。

## 可复核证据

- `delegated/paper_baselines/PAPER_CONTRACTS.json`：2016/2019 原文物理合同及缺失项。
- `delegated/paper_baselines/SCHOLARQA_COLLECT.json`：实际429搜索负回执。
- `delegated/paper_baselines/SCHOLARQA_VERIFY.json`：五个 DOI/arXiv 身份验证成功。
- `delegated/paper_baselines/ADDITIONAL_CITATION_VERIFY.json`、`LUCIE_VERIFY.json`、`CIE_NETWORK_FBE_VERIFY.json`：其余直接先例的独立身份验证。
- `delegated/paper_baselines/PRIMARY_DOWNLOADS.json`：网络下载成功/失败分别保留；作者 DLSOM PDF 已经通过网页全文读取，即使本地下载遭连接重置。
- `PAPER_BENCHMARKS.md`：本轮实际数值结果和原论文合同差异。

此报告不是模型全球新颖性认证。上表区分完整原文、摘要与本线程解释，未把引用数或题名相似视作贡献证据。
