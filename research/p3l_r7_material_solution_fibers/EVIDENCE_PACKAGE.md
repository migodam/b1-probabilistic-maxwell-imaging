# R7阅读与无损证据

[现有B1的R7 Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3lr7-material-fibers-20261028)提供完整原始证据、manifest、独立校验收据及GPT Pro阅读包。

- `Research_P3LR7_Results_and_GPT_Pro_Handoff.zip`：中文报告、源码、图、compact八对象数值与必要复核入口；不是全部过程数据。
- `R7_Complete_Evidence.tar.zst`：所有新增科学数组、封存预测、逐trial current/field/J、账本、失败与必要原始复用输入。恢复后才能按原相对路径np.load。
- `R7_Evidence_Manifest.json`：逐文件尺寸、BLAKE2b、时间与来源。
- `R7_Delivery_Verification.json`：无损逐字节恢复核验；最终匿名读回另有收据。

范围只包含本轮R7及明确必要的约10.2MB历史复用输入，不重复上传R6九GB完整包、历史训练数据或虚拟环境。43个必要复用文件复制后逐字节相等，ASSET_MANIFEST记录原始路径与恢复层次；没有保存LU矩阵，保留current/field/J/T和可重建状态。

原绝对路径是来源证据；在另一台机器恢复到新的目录层次即可，不覆盖不同字节的已有文件。纯数组审计支持新输出目录与路径重定位。归档保留原文件逻辑字节与科学字段，打包、发布和读回0新增Maxwell RHS，不提高任何科学资格。没有SHA-256附加检查，没有删除原始科学数据。
