# 单次有界开发修正

原两起点按测量有限库分别选择disk和ellipse，但后者只有pi/2方向；合法ellipse0分支被遗漏。全阶段尚余30次FE，本修正把同一公开ROI中心ellipse0起点给予known-wall、joint-cell及effective三臂，每臂最多10次实际F评价。初值不依赖目标真值，数据/目标函数/正则/材料盒/质量门不变。

属于已暴露开发适应，不作新独立确认，不修改原失败。计入原E1硬上限300，不能再扩大。触限仍NONCONVERGED。
