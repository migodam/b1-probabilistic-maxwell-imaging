# R11-S3 从这里开始

**视觉资格补充：P目标3本次离线构造丢失中横，标为 `DECODER_REPRESENTABILITY_FAIL_VISUAL_THIS_FIT`。原数值fit门仍PASS，原数据不改；不能将P称为九个清晰数字。见[补充说明](REPORT_QUALIFICATION_UPDATE.md)。**

已完成N/P各九个数字、四个物理后端共72次真实逆散射，不是模板分类。主线程已先盲读图像后解封GT。

**整轮资源协议FAIL_CAP_EXCEEDED：28,727全阶LU超过12,000。已停止新增物理，不回写预算。** 数值结果与失败仍完整保留。

1. [中文完整报告](RESEARCH_REPORT_ZH.md)，含两套九数字对照图。
2. [逐数字盲视觉审查](DIGITS_VISUAL_REVIEW_ZH.md)，特别检查8的两个孔。
3. [A精确支撑](agent_a/REPORT_ZH.md)：稀疏数字完整J核节时63.54%，不等于全成像收益。
4. [B固定总场](agent_b/REPORT_ZH.md)、[C统一损耗界](agent_c/REPORT_ZH.md)：可运行但未加速。
5. [方法CSV](METHOD_COMPARISON.csv)、[逐对象CSV](DIGIT_PER_NUMBER_RESULTS.csv)、[资源汇总](COST_SUMMARY.json)。
6. [复现说明](REPRODUCE.md)、[科学决策](GO_NO_GO.md)、[下一轮交接](GPT_PRO_HANDOFF/START_HERE.md)。

N为新258维占比模型；P为原66维decoder中高范数的设计真值。两者不能互称独立P1统计泛化。Full66 posterior与精确材料纤维均未认证。