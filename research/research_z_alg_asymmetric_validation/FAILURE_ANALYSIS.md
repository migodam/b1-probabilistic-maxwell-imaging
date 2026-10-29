# 失败、修复、边界与未执行

1. E0原物理20RHS全部成功；写入np.bool_失败。原失败job/log/源码保留；np.generic.item修复只影响序列化。通过A/RHS/解精确缓存匹配恢复完整命名空间，0新RHS，不伪称新独立solve。
2. E2已完成33启动和第一个循环点时，响应与resolvent的sigma_min键重名。保留775RHS及INTERRUPTED.json，改名response_sigma_min；复用原33启动和首节点，只付870RHS完成其余循环。阈值、种子、方法不变。
3. E1二法各10个目标未合格，包含材料越界与边界/正确器停滞；不删除失败，不把它们算成低费用成功。33多启动仅31严格合格；其余2保留。
4. 小图册区间bound在较大box太松，按事前固定halving序列得到更小认证box。不能用这些小box认证0.05积分片或径向全路径。
5. 主cost窗口看似43.3%节时，但轨迹保存不匹配。公平buffered完整trace窗口总节时仅13.6%：保留旧数值并撤回其独立总成本主张，不改20%门。记录和代数变量膨胀会抵消RHS优势。
6. 六个70dps合并snapshot原action写END序号。原文件及ledger未改，HIGH_PRECISION_SNAPSHOT_LINK_ADDENDUM与追加事件明确关联两源已支付BEGIN；无费用丢失或额外物理。
7. 补存current中间c/j/E通过纯代数重放374目标，原端点一致，0新RHS；重放CPU另记，不把它当新发现或加速。
8. 文献两查询四操作429、CSI原文captcha、绘图system Python缺NumPy等行政失败保留。文献身份核验和可得全文层级分开；绘图换已有环境，不影响物理。

NOT_RUN：完整Full66/576研究、严格Hankel或连续PDE证书、全域根数/组件覆盖、复杂monodromy、新噪声/采集泛化、完整新模型带噪后验。没有通过降低容差、改noise、换先验或删接收救结果。
