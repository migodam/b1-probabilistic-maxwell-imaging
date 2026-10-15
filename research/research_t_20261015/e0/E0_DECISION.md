# E0：物理合同可复原，优先重新取得合格基础材料

新增 full-state=0，新增 full-solver RHS=0；没有生成Green矩阵、分解或求解。采用保存的原始24×24 TM/VIE几何、8源、32接收器、complex128与原噪声。32个条件均为历史已暴露 development；其中8个有Research A复跑，其余24个明示A NOT_RUN，不记失败。

原N3在完整32条件有2/32满足Research A冻结的观测相容阈值；固定A八条件中为2/8。Research A做了一步GN后的基础图为3/8。所有8个A基础图均不是原展示N3，不能把A的残差错贴给原N3。

冻结阈值τ=23.785449004057032=sqrt(χ²_512(0.95))。σ是每complex通道RMS，实/虚部SD=σ/√2；逐源先Re后Im打包并乘√2/σ，共512个实通道；gain固定为1，各源共享同一材料。阈值来自Research A后加的重建相容规则，不是声称A30历史原本以此认证基础重建。

|case|原N3残差|A基础残差|真值噪声残差|A基础数据状态|
|---|---:|---:|---:|---|
|p250065_weak_a1_snr20|28.034654|24.198733|24.090925|INVALID|
|p250065_ring_a1_snr20|32.987001|27.065436|22.534077|INVALID|
|p250066_weak_a1_snr20|27.152058|NOT_RUN|21.913376|NOT_AUDITED|
|p250066_ring_a1_snr20|39.391876|NOT_RUN|21.570927|NOT_AUDITED|
|p250067_weak_a1_snr20|28.995031|NOT_RUN|21.766454|NOT_AUDITED|
|p250067_ring_a1_snr20|39.750525|NOT_RUN|23.835806|NOT_AUDITED|
|p250068_weak_a1_snr20|23.945286|NOT_RUN|22.542662|NOT_AUDITED|
|p250068_ring_a1_snr20|33.486828|NOT_RUN|22.956691|NOT_AUDITED|
|p250069_weak_a1_snr20|23.212852|20.929165|21.841366|VALID_NUMERICAL|
|p250069_ring_a1_snr20|31.744707|26.376896|22.434088|INVALID|
|p250070_weak_a1_snr20|33.333055|NOT_RUN|22.238503|NOT_AUDITED|
|p250070_ring_a1_snr20|67.688159|NOT_RUN|23.182112|NOT_AUDITED|
|p250071_weak_a1_snr20|31.267097|NOT_RUN|21.980030|NOT_AUDITED|
|p250071_ring_a1_snr20|65.779539|NOT_RUN|21.956017|NOT_AUDITED|
|p250072_weak_a1_snr20|38.742258|NOT_RUN|23.834282|NOT_AUDITED|
|p250072_ring_a1_snr20|49.558641|NOT_RUN|22.598620|NOT_AUDITED|
|p250073_weak_a1_snr20|23.051898|20.891919|21.718021|VALID_NUMERICAL|
|p250073_ring_a1_snr20|37.895938|25.115852|22.344127|INVALID|
|p250074_weak_a1_snr20|28.779963|NOT_RUN|23.055600|NOT_AUDITED|
|p250074_ring_a1_snr20|44.383090|NOT_RUN|22.511359|NOT_AUDITED|
|p250075_weak_a1_snr20|30.296817|NOT_RUN|22.735577|NOT_AUDITED|
|p250075_ring_a1_snr20|40.899869|NOT_RUN|23.171442|NOT_AUDITED|
|p250076_weak_a1_snr20|29.615177|NOT_RUN|23.673932|NOT_AUDITED|
|p250076_ring_a1_snr20|35.662134|NOT_RUN|23.149688|NOT_AUDITED|
|p250077_weak_a1_snr20|24.640703|22.994968|23.047116|VALID_NUMERICAL|
|p250077_ring_a1_snr20|31.544581|26.334774|22.654196|INVALID|
|p250078_weak_a1_snr20|25.973154|NOT_RUN|21.983053|NOT_AUDITED|
|p250078_ring_a1_snr20|41.455241|NOT_RUN|23.235433|NOT_AUDITED|
|p250079_weak_a1_snr20|27.708472|NOT_RUN|23.847844|NOT_AUDITED|
|p250079_ring_a1_snr20|35.797072|NOT_RUN|21.232579|NOT_AUDITED|
|p250080_weak_a1_snr20|30.738861|NOT_RUN|23.002724|NOT_AUDITED|
|p250080_ring_a1_snr20|42.869752|NOT_RUN|21.130811|NOT_AUDITED|

原验证档的clean response与A28对应真值full prediction逐位一致；原S作用于保存current同样一致。保存的场方程、材料本构与预测检查最大相对误差1.97e-15，没有已有模型不匹配证据。元数据噪声种子重新生成的观测误差为0。

原N3 pool材料与已付费cache材料不是逐位相等：mass打包往返带来至多4.45e-16的ULP量级差异，各case另记录exact=False和numerical_match=True。这与Research A一步GN造成0.028–0.068量级材料变化是两件事，不能混淆。缓存的matching endpoint按既有数值容差解释，不声称byte identity。

首个weak条件真实噪声白化范数24.090925超过95%阈值，属于原噪声模型允许的小概率实现。不得因此宣称模型错误，不得为取得通过而改σ/阈值。即使真实材料恰好不满足这一噪声球，也不排除其他拟合噪声的合法材料存在。

E0没有阻断本轮已授权的E0-Q。应使用原数据/材料盒作更充分的全材料修正，明确保存KKT、线搜索、active bounds、每个匹配full endpoint和成本。不能把当前未合格基础图交给G冒充合格初始逆解。

原始N3缓存有完整fields/current并已零求解检查；A候选只保存material/full prediction，细节见missing manifest。`VALID_NUMERICAL`不含严格输出误差界。
