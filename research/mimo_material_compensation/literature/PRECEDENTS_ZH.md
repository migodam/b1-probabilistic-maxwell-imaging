# 本轮主张与直接先例

| 原始工作 | 已给出的结论 | 本轮仍须自己证明的事情 |
|---|---|---|
| Bucci、Franceschetti；Bucci、Isernia | 散射场准带限，系统孔径与角度采样限制可读自由度 | 对有限面积任务，合法材料能沿同一带噪数据球变化多远。由最坏端口功率、标定先验和完整材料对给次频读数的充分性。 [1987 原文](https://web.eng.ucsd.edu/~massimo/ECE287C/Handouts_files/Bucci-Franceshetti_1987.pdf) · [1997 Radio Science](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/97RS01826) |
| Bonnet-Ben Dhia、Chesnel、Nazarov | 已在自由空间有限入射与接收方向下构造不可见的穿透介质，定义相对散射矩阵及其核 | 在有损端口、实非负介质、指定 ROI 面积和幅度下，严格给出合成读数、补偿预算和可检测边界。 [原文](https://arxiv.org/abs/1410.8382) |
| Bera、Bonnet-Ben Dhia、Chesnel | 波导中用材料导数延拓完全或相对不可见解，单模有完整证明，多模提供计算构造 | 对给定带噪读数给出原生二维网格的整段严格宽度，并区分精确等场与噪声球内的连续路径。 [原文](https://arxiv.org/abs/2005.06939) |
| Ammari、Kang、Lee、Lim及合作者 | 层状材料可使若干散射系数消失、增强近隐身；电磁矢量情形亦有先例 | 同向调节两件非负材料的外区、在指定强散射加载阵列下审计一个被动探头的全部负载和所有波束。构造环层本身已有明确前例。 [Helmholtz](https://arxiv.org/abs/1110.1922) · [Maxwell](https://arxiv.org/abs/1212.5685) |
| Bouchet 等 | 给定两个散射配置和相干激励，差分散射矩阵的首奇异方向是二元辨别强基线，有实验实证 | 对指定材料类的最坏合法外区再拟合。受限在同一硬件时，比较任意相干波束的必要能量与另一频率整类测试的充分能量。 [PRL 原文](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.127.253902) |
| del Hougne | 通过已知可调负载读取不可直接接入的端口；三个负载、互易符号歧义和解除歧义的额外耦合负载已有理论和实测 | 在新增探头原本借以分辨材料的情形，求其全部稳态复负载面对未知材料补偿的能力上限。三态反演隐藏端口参数本身无原创主张。 [Virtual VNA](https://arxiv.org/abs/2403.08074) · [Virtual VNA 2.0](https://arxiv.org/abs/2409.10977) |
| Atkinson、Fedorov | 可用竞争模型最小拟合误差选下一次实验，涵盖非线性回归和序贯设计 | 对有界强多次散射模型给严格方向余项，并纳入会重新补偿的两份材料。留一拓扑结果没有显示新指标普遍优于这一类判别思想。 [Biometrika 原文](https://doi.org/10.1093/biomet/62.1.57) |
| Alberti、Santacesaria | 合适的有限维系数先验与边界测量可建立唯一性、稳定性、重建 | 此处硬件、空间激励、噪声、能量和冲突两类都事先给定；宽先验中报告实际冲突材料。 [Forum of Mathematics, Sigma 原文](https://arxiv.org/abs/1803.04224) |
| Griesmaier、Harrach | 在亥姆霍兹散射的远场算子上建立模有限维单调性；更正收紧了形状结果的假设 | 有限电路端口和接收噪声下的任务量与覆盖宽度需另行验证。 [原论文](https://epubs.siam.org/doi/10.1137/18M1171679) · [作者更正](https://www.math.uni-frankfurt.de/~harrach/publications/farfield_monotonicity_erratum.pdf) |
| Chen；Zhong、Chen；Zhong 等 | 从电流空间分解或重写积分方程改善求解。原模型物理与充分迭代的 CIE 是重要物理对照 | 本轮在测量空间除去可补偿响应，并以原方程验证多份物理解。将单个求解器的失败误作排除会直接给出错误的材料结论。 [SOM](https://doi.org/10.1109/TGRS.2009.2025122) · [TSOM](https://doi.org/10.1088/0266-5611/25/8/085003) · [NIE](https://doi.org/10.1109/TAP.2016.2535492) · [CIE](https://www.mdpi.com/2313-433X/5/2/27) |
| Boyd、El Ghaoui、Feron、Balakrishnan | 线性系统的有界实、耗散与矩阵不等式已有完备工具 | 原生全波材料对导出的两个反馈状态、全负载上界、另一频率的整类下界以及能量量词须共同核对。 [SIAM 专著作者页](https://web.stanford.edu/~boyd/lmibook/) |

上述先例均有原论文或作者平台可核对。材料路径证明本身不建立“首次构造不可见材料”；单对最优波前、负载消元、线性盒支持和有界实引理各有成熟来源。可争取的贡献集中在一个固定硬件上，指定面积、合法材料预算、噪声与测量协议的同时定量上下界。
