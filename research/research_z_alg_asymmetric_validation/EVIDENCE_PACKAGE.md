# 阅读与完整实验原始证据

本轮公开在既有B1，不新建仓库。Release目的地为 [Z-ALG非对称验证](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/z-alg-asymmetric-validation-20261029)。实际发布/匿名回读另有独立收据，不能由这个预声明链接推断发布已完成。

- Research_Z_ALG_E1_Complete_Evidence.zip：全部本轮物理/代数实验数组、完整轨迹、原始账本、失败、配置、源码、原用户理论输入、图与独立审阅。
- Research_Z_ALG_E1_GPT_Pro_Handoff.zip：精简中文报告、图、源码、主要compact结果与Pro接口；不是全部逐trial过程。
- Z_ALG_Evidence_Manifest.json及Z_ALG_Delivery_Verification.json：逐文件逻辑字节、BLAKE2b、时间、ZIP每条目逐字节回读及排除项。

普通NPZ可直接np.load(allow_pickle=False)，没有dtype下转、科学字段删除或有损压缩。所有非合格结果和两次报告错误保留。

外部检索取得的第三方论文全文/网页/长摘要只在本地保存，公开包改用原出版社链接、元数据与自行撰写的简明证据账本；它们不是本轮新实验原始数据。此排除不删除本地文件，也不影响全部新实验数组与失败证据。原用户提供的Z-ALG报告/源码包作为输入保留，旧R6/R7大证据不重复打包。
