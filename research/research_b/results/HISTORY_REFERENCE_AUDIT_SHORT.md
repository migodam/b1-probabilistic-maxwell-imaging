# 固定 joint-history 数值审计

4个固定parent、各2noise；32政策标签行，16个唯一history。原proper coefficient prior + persistent bounded gain；只用已选观测。

- finite2048 对17阶最大event质量差：0.186853037
- finite2048 对17阶最大mean/prior-width差：0.121553777
- finite2048 对17阶最大已选candidate expected-L1-risk差：0.0162756365
- 13/17最大event差：0.0907806963
- 13/17最大mean/prior-width差：0.11889838
- 13/17最大SD相对差：3.04713586
- UNRESOLVED唯一history：12 / 16
- 峰值RSS：1216.1 MiB；新增Maxwell0、GPU0。

这仅检查固定子集、同一连续proxy target的数值近似。未增加21阶/CDF/全案例campaign；17阶本身仍只由13/17比较支持。
原算法、Selected-MAE统计和候选决策保持原样；参考event mass、expected risk与校准主张必须同时披露上述误差及未稳项。主会话保留最终科研裁决。
