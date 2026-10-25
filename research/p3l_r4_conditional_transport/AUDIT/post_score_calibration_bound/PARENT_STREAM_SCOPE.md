# 封账后源代码核对：抽样流与同时下界

实际冻结evaluator并非每份parent反复使用同一个seed。`tools/evaluate_r4_density.py:154–161`用`SeedSequence([sample_seed,index])`生成parent级seed，再传给full66 predictor；`src/p3r4/density.py:446–468`使用该CPU Generator生成K128 base draws。对保存的两份FROZEN_DENSITY metadata只做整数种子检查：512/256各自均无重复seed。未重新推断、读取GT数组或求解Maxwell。

方法之间复用同一parent随机流是common-random-number控制；两批也可能复用index流。Bonferroni同时覆盖不要求52项检验彼此独立，只要求每项二项区间在其批内成立。正式TV/KL解释仍依赖理想IID随机流、固定模型与原离散noise合同，不能将无重复整数seed当作数学独立性证明，也不能宣称条件于冻结PRNG比特的位级置信证书。原始分数、门槛与独立审阅没有改写。

该核对晚于阅读ZIP/后补ZIP生成；它验证已有代码的抽样身份，后来的出版提交单独保存，不伪称原包中已有。
