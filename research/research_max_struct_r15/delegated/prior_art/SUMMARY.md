# R15 closest-prior：有限证据摘要

本轮保留六篇直接相关原始论文，未做系统综述、核心证明或R15新颖性/有效性判决。检索与核验日期为2026-11-16（Asia/Singapore）；无physics、GT、远端任务或环境安装。

早期SOM的原论文已把signal-subspace取得的信息与非线性优化的余部结合；Twofold2009摘要进一步明确利用内部current-to-field谱，先在低维current空间产生初值，再扩大current空间。这里证据只到摘要，不能据此声称已逐式核对固定19维PGP/GP材料分母或新的增量更新。[Chen2009](https://doi.org/10.1364/JOSAA.26.001022)，[Zhong–Chen2009](https://doi.org/10.1088/0266-5611/25/8/085003)。

两投影的四个交集通常不构成全部空间：Halmos在开头另列generic remainder，Theorem2在此部分保留CS耦合块。Björck–Golub则提供principal angles/vectors的QR/SVD计算框架。因此文献对照应保留generic/主角度信息，不能把“四交集标签”当作一般非交换投影的完整分解。这是标准子空间几何，不是材料可辨识性结论。[Halmos1969](https://doi.org/10.1090/S0002-9947-1969-0251519-5)，[Björck–Golub1973](https://doi.org/10.1090/S0025-5718-1973-0348991-3)。

部分/块DR已有原始理论，但Combettes–Pesquet要求具体算子类别、解存在、抽样/误差及正激活等条件。Elser的RRR完整更新具有full残差零与solution/shadow的关系，局部flow与搜索阶段分开讨论。将这些结论用于“只投影increment”的具体方案，仍须核实更新映射及固定点关系，不能只凭DR/RRR名称移植；本记录不据此判定候选方案成立或失败。[Combettes–Pesquet2015](https://doi.org/10.1137/140971233)，[Elser2021](https://doi.org/10.1186/s13663-021-00697-1)。

## 已核对的六条引用

|键|原始论文及作者|出版记录|DOI|
|---|---|---|---|
|P1|Xudong Chen, Application of signal-subspace and optimization methods in reconstructing extended scatterers|JOSA A26(4),1022–1026,2009|[10.1364/JOSAA.26.001022](https://doi.org/10.1364/JOSAA.26.001022)|
|P2|Yu Zhong; Xudong Chen, Twofold subspace-based optimization method for solving inverse scattering problems|Inverse Problems25(8),085003,2009|[10.1088/0266-5611/25/8/085003](https://doi.org/10.1088/0266-5611/25/8/085003)|
|P3|P. R. Halmos, Two subspaces|Trans. AMS144,381–389,1969|[10.1090/S0002-9947-1969-0251519-5](https://doi.org/10.1090/S0002-9947-1969-0251519-5)|
|P4|Åke Björck; Gene H. Golub, Numerical methods for computing angles between linear subspaces|Math. Comput.27(123),579–594,1973|[10.1090/S0025-5718-1973-0348991-3](https://doi.org/10.1090/S0025-5718-1973-0348991-3)|
|P5|Patrick L. Combettes; Jean-Christophe Pesquet, Stochastic Quasi-Fejér Block-Coordinate Fixed Point Iterations with Random Sweeping|SIAM J. Optim.25(2),1221–1248,2015|[10.1137/140971233](https://doi.org/10.1137/140971233)|
|P6|Veit Elser, Learning without loss|Fixed Point Theory Algorithms Sci. Eng.2021,12|[10.1186/s13663-021-00697-1](https://doi.org/10.1186/s13663-021-00697-1)|

## 证据边界与检索留账

- ScholarQA CLI一次3-query小bundle：六个paper/snippet请求中五个429；未重试扩搜。一次最终verify解析5/6 DOI，Halmos的AMS DOI在S2未解析；其身份由AMS原文索引及出版社沉积的Crossref记录另核，未编造S2 ID。
- Undermind只做定向metadata lookup/get-info，未启动deep search。其Halmos2010、Björck–Golub1971及两个arXiv首稿年份不能当作期刊出版年份；按原出版社/DOI记录分别记1969、1973、2015、2021。P5预印本首稿2014，P6首稿2019，均与期刊年分开。
- P1/P2只有原出版方摘要（Tier B）。P3/P4读取到可定位的AMS原文索引段落/摘要，直接PDF403，因此标partial，不声称完整审读。P5原作者[arXiv v4](https://arxiv.org/html/1404.7536v4)及P6[出版社开放全文](https://link.springer.com/article/10.1186/s13663-021-00697-1)已核对上述指定位置。
- 引用入口复用R10/R12/R13已有closest-prior文档；这些内部文档仅作检索线索，实质主张以这里列出的primary论文为证。不绕过付费墙、不安装工具；未取得全文不当作没有前人覆盖的证据。

逐条支持、限定和推断见`CLAIM_LEDGER.csv`，检索失败与元数据冲突见`RETRIEVAL_RECEIPT.json`。具体空间/更新规则的数学有效性、是否构成增量及最后科学结论留给父线程。
