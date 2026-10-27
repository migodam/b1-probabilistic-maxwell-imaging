# 阅读包构建状态

当前仅完成构建器和`READING_ZIP_PLAN.json`清单，**尚未生成ZIP**。根线程在最终科学提交之后，使用该提交的完整ID显式执行：

```sh
../../.venv_nn/bin/python research/p3l_r6_task_optimal_maxwell/GPT_PRO_HANDOFF/build_reading_zip.py --build --source-commit <final-source-commit>
```

命令在R6仓库根执行。输出固定为Inv_SLAM根下`research/delegated/r6_delivery_20261027/Research_P3LR6_Results_and_GPT_Pro_Handoff.zip`。构建器要求所给提交等于当前HEAD，关键最终报告和构建器已经提交且没有未提交修改；已存在ZIP不会被覆盖。没有`--build`时只刷新清单。构建后逐文件解包做字节核对，并写`READING_ZIP_RECEIPT.json`。

包内保留R6相对目录结构，收录最终中文报告、理论、阶段报告、小JSON、关键CSV、图、实际Python来源、静态网页，以及主D两份无损数值NPZ和映射清单。每条件的全部256对象、六动作、三seed均在紧凑包中。清单排除大raw、输入历史ZIP、J/H与field/current、优化trace、训练模型、逐对象context和虚拟环境；原始证据通过`RAW_RELEASE_REFERENCE`里的实际稳定分包manifest/receipt引用。发布模板不等于已发布URL，最终在线入口以根线程发布回执为准。

这是给下一轮GPT Pro的阅读包，不能代替完整原始证据恢复包，也不能将已经暴露的GT再称为未见评测。条件分开，点动作不是posterior样本；最终科学判定以根线程的报告和GO/NO-GO为准。本构建器不会训练、推断、评分或调用物理。
