# 源码与完整无损证据

本轮科学实验和本地交付已封账。公开发布及匿名核验另记，不以文件压缩成功代替科学或发布资格。

四个tar.zst分卷共 **1,714,274,978 bytes（约1.714GB）**，均小于2GB。全部121,850个raw文件包含在分卷中，16个冻结checkpoint及所需旧源码/背景阵列单独一包；原始文件保留，没有删除或降精度。每个tar成员已解压并与原文件逐字节比较通过。

|分卷|文件数|压缩bytes|本地无损核验|
|---|---:|---:|---|
|[test-a-raw.tar.zst](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/p3lr5-b-confirm-20261027/test-a-raw.tar.zst)|59,119|798,479,866|PASS|
|[test-b-raw.tar.zst](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/p3lr5-b-confirm-20261027/test-b-raw.tar.zst)|60,607|799,818,128|PASS|
|[shared-probe-scores-ledger.tar.zst](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/p3lr5-b-confirm-20261027/shared-probe-scores-ledger.tar.zst)|2,126|103,673,976|PASS|
|[essential-reused-source-checkpoints.tar.zst](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/p3lr5-b-confirm-20261027/essential-reused-source-checkpoints.tar.zst)|253|12,303,008|PASS|

对应同名JSON manifest随Release提供；记录相对于Gaussian目录的安全相对路径、大小和来源。原checkpoint另有BLAKE2b内容清单；没有SHA256检查。归档原样保留complex128/float64与网络原float32输出，不是有损量化。恢复到全新隔离目录，不能直接覆盖历史工作区。实际恢复路径以manifest为准。

[完整Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3lr5-b-confirm-20261027) · [Pro阅读ZIP](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/p3lr5-b-confirm-20261027/Research_P3LR5_B_Confirm_Results_and_Pro_Handoff.zip) · [复现说明](REPRODUCE.md)

阅读ZIP只含报告、代码、冻结配置、紧凑指标与真实图，不冒充完整原始证据。新原始GT、材料/预测、完整RHS账本、失败日志、端点与轨迹在四分卷中。旧32768 TRAIN/2048 VAL不重复打包，历史训练的完整来源仍见原R5公开Release；本次预测不需要重新训练。

实际包操作wall合计101.27秒、CPU40.92秒；中途一次“归档已验证、后续记账字段重复”的异常保留，修复并复用已验证A包，没有重算物理或改动数组。次级交付总wall和未知小额启动费用见LOCAL_EVIDENCE_PACK_RECEIPT，不虚称完整CPU全归因。所有交付成本不计为Maxwell RHS。
