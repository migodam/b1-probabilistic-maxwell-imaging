# 网站表格与图集来源

本次仅呈现已有结果，无新评分、优化或科学判断。

| 展示 | 实际来源 | 精确行选择 |
|---|---|---|
| 六方法48对象主表 | `analysis/stack_comparison/METHOD_VIEW_SUMMARY.csv` | cohort=uniform48，initialization=uniform，view=primary；每solver唯一一行，available_cases=NRMSE_n=48 |
| 六方法共同CPU列 | `COMMON_ACTUAL_CPU_SIX_AND_TWELVE_SUMMARY.csv`（同目录） | uniform48 / uniform / common_actual_CPU_six / uniform_six |
| NN16前两列 | METHOD_VIEW_SUMMARY.csv | factorial16 / uniform或NN / primary；每组合唯一一行，available_cases=NRMSE_n=16 |
| NN16共同CPU两列 | 共同CPU summary CSV | factorial16 / uniform或NN / common_actual_CPU_twelve / factorial_twelve；不混用独立six预算 |
| 主表实际均值CPU | `RESEARCH_REPORT_ZH.md` 第一张六方法表 | 原正文字符串，未重新平均 |
| 中文HTML报告 | `RESEARCH_REPORT_ZH.md` | 全文转换，TeX公式保留可读源码；未安装依赖或引入外部脚本 |

网站数值按CSV原值格式化为五位小数；机器可读的精确原行保存在 `SITE_TABLE_SOURCE_AUDIT.json`。六solver顺序固定为CIE、Dual-CIE、G19、DR、A+-MM/R13、A+-Arc/R7-readout，不新增排名或阈值。

完整图94张：五方法绝对色标12、逐行诊断12、五方法曲线2、Arc/MM全48共8、NN全16共8、NN曲线2、动力学2、native48。图片与PDF只链接实际存在路径，原文件名及图中完整ID保留；native按完整对象ID展开。没有private_operational或fixture链接。

运行 `python src/build_supplement_site.py` 生成入口、HTML报告及8组图集；收据见 `SITE_BUILD_RECEIPT.json`。重建时只收录已经存在的文档与实际 `DOWNLOADS.json`；不存在清单时不生成猜测的下载链接。此次未发布、未改原R13首页或科学代码；尚未进行浏览器视觉检查。
