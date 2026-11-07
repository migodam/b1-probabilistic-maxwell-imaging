# Green-current / ordinary-material曝光DEV因子消融

按预声明只在80确认结束后执行3个公开DEV，原始prior和has_reg完整保留，普通材料直接进入继承的global_dual `_baseline_readout` 分支。已验证H本身未改，无额外恒等式测试。3预测先封存后读旧曝光GT。

|case|scalar+ordinary NRMSE|scalar+Green NRMSE|Green+ordinary NRMSE|Green+Green NRMSE|
|---|---:|---:|---:|---:|
|r10_twod_0001|0.422556|0.098895|0.422033|0.105166|
|r10_twod_0004|0.110322|0.515137|0.515136|0.515137|
|r10_twod_0022|0.096292|0.337003|0.554621|0.100878|

Green-current + ordinary-material在三例均材料FAIL且CN99 FAIL；桥corrector到3000上限。这个结果不支持把BothBlocks曝光桥成功完全归给current单独替换。完整轨迹和选点也发生改变，因此不是独立的同状态因果证明；新16确认结果仍必须单独保留。

3个新DEV child的精确wait4 CPU总和128.297977秒，失败0、重试0，成本不合并入80确认。源/配置启动前快照、完整receipt、封存、scores及cost ledger均在本目录。
