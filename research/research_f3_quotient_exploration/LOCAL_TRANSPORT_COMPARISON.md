# 四臂的实际有限传输比较

两条件均通过新stationarity检查，因此A/B/C/D均实际运行。均在相同原材料支持、noise/gain、源/接收器和初始材料方向下，申请s与2s两个材料距离。

A冻结商坐标切向；B每步更新prior tangent和Jacobian并对齐；C对材料目标做合法retraction；D加入prior法向曲率、一次真实Maxwell方向Hessian，并在第二目标使用第一次full响应误差及刷新Jacobian。D是一个冻结的小型二阶/测量反馈启发式，不是完整非线性全局corrector，也不带最优保证。

| 条件 | 方法 | 目标RMS | 实际RMS | ΔΦ | 密度正值 | 距离匹配 |
|---|---|---:|---:|---:|---|---|
| ring | A_frozen | 0.0052756 | 0.0052756 | 20.633656 | True | True |
| ring | A_frozen | 0.0105511 | 0.0101229 | 80.586403 | True | False |
| ring | B_updated | 0.0052756 | 0.0052756 | 20.633656 | True | True |
| ring | B_updated | 0.0105511 | 0.0054075 | 21.450607 | True | False |
| ring | C_native | 0.0052756 | 0.0052756 | 5.196378 | False | True |
| ring | C_native | 0.0105511 | 0.0105511 | 17.032442 | False | True |
| ring | D_intrinsic_feedback | 0.0052756 | 0.0052756 | 4.502304 | True | True |
| ring | D_intrinsic_feedback | 0.0105511 | 0.0105511 | 36.853644 | False | True |
| weak | A_frozen | 0.0052756 | 0.0052756 | 10.643981 | True | True |
| weak | A_frozen | 0.0105511 | 0.0105511 | 51.105662 | True | True |
| weak | B_updated | 0.0052756 | 0.0052756 | 10.643981 | True | True |
| weak | B_updated | 0.0105511 | 0.0065312 | 15.645211 | True | False |
| weak | C_native | 0.0052756 | 0.0052756 | 1.140604 | True | True |
| weak | C_native | 0.0105511 | 0.0105511 | 3.098031 | True | True |
| weak | D_intrinsic_feedback | 0.0052756 | 0.0052756 | 1.632730 | True | True |
| weak | D_intrinsic_feedback | 0.0105511 | 0.0105511 | 4.230149 | True | True |


primary s处四臂均匹配。weak的ΔΦ从冻结10.644降到C的1.1406，D为1.6327；ring从20.634降到C的5.1964，D为4.5023。材料坐标/retraction有有限步收益；高阶D的改进不稳定，并未增加任务尺度competitive成功数。

B的第二目标在两条件都未匹配，ring A的第二目标也未匹配。它们不能被当作同2s的有效比较或免费成功。D/C第二目标匹配，但weak ΔΦ3.098/4.230都不competitive。

D使用一次8源方向solve求Maxwell二阶项；旧缺失Hessian未重算。新CM/II/方向/normal correction及预测误差全部保存在新anchor的INTRINSIC与TRANSPORT_PAIR数组。两个端点不等于已知完整路径；后续中点诊断不支持所选参数直线的competitive保持性。

没有rank-deficient exact tangent证据；这里在有限噪声材料区域探索。NO-GO指当前D增量/学习主张，不是说一般Maxwell有限反馈没有价值。
