# 额外零物理数学单元

`EXTRA_MATH_UNITS.py`已运行，四组全部PASS。只导入NumPy/SciPy与复制的静态prior，新增Maxwell states/RHS/训练=0。

|单元|实测结果|科学范围|
|---|---|---|
|真实decoder三维Q×D同材料gauge|pure Q=1,D=(1.3,.105)；transition Q=.95,D=(1.368421,.110526)，材料最大差4.48e−16|3D固定子模型；不证明Full66分层正质量overlap|
|遗漏质量恒等式与任务振幅界|100个有限、已知不交partition；恒等式最大误差4.44e−16；bound最小slack .009806|有符号任务；不是未访问Maxwell域的认证|
|强方向补偿f=z+h³|sigma=.1/.05/.02均得mean(h)=.5,var=1/12,P(h>.5)=.5|z∈[-2,2]可完全补偿h∈[0,1]；非零c不推出新增收缩|
|Gaussian derivative covariance identity|含atom的cubic toy，均值/二阶矩/事件三任务；中心差分step1e−5最大误差3.65e−10|有限求积Q1；无定向舍入，无物理结论|

附错误q反例：y无关的三点q其真实teacher derivative为0，而Cov_q(M,F)/sigma²=4.1667。任意近似q的covariance不能自动作为真实teacher导数。

所有结果与参数见 `EXTRA_MATH_RESULTS.json`。这补充所列四类单元，不将其余未运行的P0区间/3D密度积分检查自动改PASS。
