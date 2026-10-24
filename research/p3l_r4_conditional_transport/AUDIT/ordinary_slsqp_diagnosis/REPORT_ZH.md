# 原ordinary SLSQP：已付费试探诊断

{
  "status": "SAVED_ONLY_ORDINARY_DIAGNOSIS_COMPLETE",
  "parents": 64,
  "paid_trials": 1024,
  "all_fallback": 64,
  "zero_callback_parents": 0,
  "total_callbacks": 198,
  "feasible_callbacks": 0,
  "selected_callbacks": 0,
  "callback_Phi_improvements": 0,
  "first_trial_step_norm": {
    "minimum": 0.9748915427096035,
    "median": 2.690297540323334,
    "maximum": 9.585083303265597
  },
  "first_trial_direction_cosine": {
    "minimum": 0.9999999999999998,
    "median": 1.0,
    "maximum": 1.0000000000000004
  },
  "first_trial_Phi": {
    "minimum": 290.9032693313363,
    "median": 4632.765658292831,
    "maximum": 179106.00998204813
  },
  "first_minus_gradient_error": {
    "minimum": 2.436383457711112e-16,
    "median": 2.3592785926533365e-14,
    "maximum": 1.4223628423238166e-13
  },
  "saved_anchor_gradient": "Jq*probitPDF then Jx.Tsqrt2/sigma*pack(F0-y),allalready-paidonline",
  "new_RHS": 0,
  "new_LU": 0,
  "new_states": 0,
  "GT_read": false,
  "mu_read": false,
  "old_run_modified": false,
  "interpretation": "SLSQPcallbackrecordsarenotproofthatobjectiveimproved; budgetstop+fallbackcannotestablishnostrengtheningpossible"
}

只读缓存在线J/原残差、原试探预测与callback事件；不读取GT或mu，不作新求解。原320case与封账不改。callback标签沿用原ACCEPTED_ITERATE，不能把有callback等同于下降或最终可用。

## callback和已验证查询的区别

安装SciPy1.18.1的包装器仅按major iteration计数递增触发callback，没有检验下降或数据可行性。因此原ACCEPTED_ITERATE名字不是内部线搜索接受证书。实际64对象有198个callback，不是零callback；所有callback的Phi都超原噪声球且不比锚点低，原callback-only选择保留64锚点。但所有64对象都有原先已经付费的有效Full查询Phi低于返回值，差值最小.00015779、中位.00207002、最大.10562763。这证明不能把返回停滞解释成数据Phi不能下降。

原始初次位移确实是负数据梯度（最大向量误差1.42e-13），范数中位2.690，最大9.585；第一Phi中位4632.77、最大179106。方向符号与缓存J链式法则相符，但单位初始步在强非线性/各向异性白化目标上很差。这里只据保存序列诊断，不把callback具体低层线搜索时序说成已证明。

另行预注册的best-validated-query补充仅在旧DEV0..15查询池执行：48case中ordinary16/16得到更低Phi，mu32case未找到比原输出更低的可行材料目标。梯度失败不使已经通过原Full/先验合法性/后向残差检查的世界失效；optimizer未接受的查询也可在明确称BEST_VALIDATED_QUERY_CANDIDATE时作为新附加输出。原callback flags、320个freshL1输出、封账和所有已付查询成本都不改，没有L1新选择器执行、GT读取或追认primary gain。

Full-J TR新DEV16两目标控制协议另冻于configs/ORDINARY_FULLJ_TR_DEV16_PROTOCOL.json，初始半径.5/上限2、ratio固定规则和16状态/256RHS明确；3项纯代数KKT/白化probit测试通过。仅在父PPC完成明确释放后运行，不自动启动新L1控制。
