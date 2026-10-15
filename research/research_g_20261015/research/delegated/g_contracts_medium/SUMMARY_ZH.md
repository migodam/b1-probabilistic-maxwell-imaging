# G 软件合同交付

18项通过：预算9项、侧车9项。仅临时账本/人工数组fixture；没有运行物理runner、solver、优化、GPU或SHA检查。真实新增states/RHS均0。未改g_budget.py（由root修补）。

API：`code/g_contract.py` 的 `validate_episode(dict, asset_root=None, check_assets=False) -> True`，违规抛 `ContractError`。`check_assets=True`只检查受控根目录中的文件和NPZ array_key，不能证明内容物理有效。

侧车：`GSearchEpisode.schema.json`；最小软件例子 `examples/GSearchEpisode.minimal.json`。版本GSearchEpisode/v1，target.id为LEGACY_FROZEN_A30_24x24_scalar_TM，packing_version为A30_SOURCE_MAJOR_RE32_IM32/v1，gain_version为FIXED_ONE/v1、gain数值1。这是A30原source-major Re32/Im32打包的包装标签，不是径向P2约定。

所有candidate路径均支持相对字符串或{path,array_key}；validation.full_validation_receipt为相对路径或资产引用；solver_diagnostics、constitutive_diagnostics及full_residual均为{value,reason}，数值时reason=null，缺失时value=null并给非空理由。schema不能单独表达observations中ID唯一和全树在线禁用键，必须同时调用API。

probability三项mass必须null；可选proposal_density={value,interpretation:PROPOSAL_ONLY}。拒绝重复observation_id、错误版本和物理身份、在线truth/oracle/hole_count/topologylabel键、绝对/逃逸路径、无理由缺失诊断及NaN/bool伪数值。

测试文件：tests/test_g_budget.py、tests/test_g_contract.py。独立输出和资源凭据分别为本目录下test_g_budget.txt/test_g_budget_RECEIPT.json、test_g_contract.txt/test_g_contract_RECEIPT.json。预算wall约0.130s/CPU0.106s；侧车wall约0.106s/CPU0.097s；BLAS线程1。

边界：当前venv无jsonschema依赖；未安装新包。schema的JSON结构/required/version与例子对应已测，实际拒绝行为通过轻量API完成测试，没有声称运行第三方Draft2020-12 validator。
