# E0 裁决记录：READY_FOR_BOUNDED_E1

E0零新增物理审计已完成，原prior/physics/noise可复原。原TRAIN生成器512条件精确重放，新增256个独立material-only parent，无新paired训练数据。历史T/G候选的noise-ball成功没有建立原P0支持；现存近点似然对中没有两端PHASE合格的pair。

主缓存screen的8064比较未给出结构合格且d_pred<=2的P0_EXACT pair；最接近结构合格距离44.8267177317206。缓存筛选不足以判断全局唯一性或排除未搜索pair。

因此可进入既定上限内E1：最多64个新full-state attempts、2048 complex RHS；所有失败/拒绝/line-search/验证/benchmark construction均收费，不借用历史T/G余额。E1只寻找并复核原P0竞争材料见证，32固定noise seeds重放与局部扰动按冻结协议执行；在线端不能读取隐藏端点、clean means或生成参数。

E0没有授权E2。只有E1至少两组不共享material parent、结构分离、重复噪声下竞争、无明显先验尾部或数值问题的P0_EXACT pair，且主线程科学复核通过，才可进入条件性E2。否则NO_GO_FOR_GENERATOR_COMPARISON，歧义仍未解决；不加训练/预算。COMMON_CORRECTOR_SPEC只是待激活合同，IMPLEMENTATION_NOT_ACTIVATED。最终E1/E2结论由主线程另行复核。
