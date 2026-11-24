# B3 合成确认数据接口

实际生成了 12 件独立介质和 48 份拟合视图。每件在 400 MHz、两个公开发收几何和两个配对噪声强度下重复，统计单位仍是独立材料。已有六个 SRAP 程序化几何家族生成八个新实例。其余四件来自两个原创建模家族，每个家族两件。生成时置换内部家族顺序，公开接口只返回 `b3m_0001` 到 `b3m_0012` 等匿名标识。对每件所有正材料格使用同一有限范围 \([2.5,4]\)，背景保持零。生成固定了几何规则和种子，适用于预先冻结的算法。 [生成源码](generation.py)

## 公开接口

`holdout.data.planned_cases()` 可以在没有实际数据的情况下枚举标识。例 `b3m_0001__pw180_10x20__n02` 表示第一件匿名资料、180° 几何、2% 噪声。`get_case(id, root=...)` 与 `load_public(oid, aperture, noise, root=...)` 返回 `(case_dict, public_arrays)`，其中 `E/B/E0,S,Y,centers,radius,k,delta,spectrum` 为求解和证书所需字段，默认逆网格 32²。`get_native_case(id)` 绑定相同读数和归档 96² 几何。`get_native_model(aperture)` 可先读几何，给同一阵列上的多份证书复用。两个几何 ID 为 `PW360_10x20` 与 `PW180_10x20`，噪声输入为 `0.02` 与 `0.10`。 [公开加载器](data.py)

拟合端口沿用原 SRAP 10 发 20 收的阵列和原 96² 二进制 FFT 算子。隐藏发射角为 \(9^\circ+36^\circ k\) 的十个角，隐藏接收角为 \(4.5^\circ+18^\circ k\) 的二十个角，分别避开两个拟合阵列的所有发射角和接收角。独立预测同时面对新源和新接收角。半环和全环各自生成一个归一化的复高斯噪声方向，2% 和 10% 共用该方向。公开半径为原名义噪声量加 \(10^{-6}\lVert Y\rVert_F\) 并向上取整。 [历史阵列来源](../../srap_cie_array_study_20261121/data2d.py)

## 原生审计和离线评分

`holdout.physics.NativeForward(spectrum, 96)` 使用原 FFT 嵌入和接收几何。`nearest_lift(m32)` 将每个 32² 值完整铺到对应 3×3 原生格。凡原生剪枝或原生残差结论，以同一延拓和 96² 的残差合同为准。32² 求解耗时和拟合值独立记录。两套离散模型的差异由原生验证体现。

`verify_truth.py` 独立用 96² 有源格 LU、经验证逆范数、区间原生 Green 系数和伴随误差界，证明保存的参考材料位于各公开误差球中。报告实际生成求解、重建几何、证明费用和任何失败。该程序只在冻结后读取离线材料。

`holdout.scoring.score(observation_id, m32, freeze_record=...)` 独立以原生 96² 正演验证拟合读数和隐藏发收对。也接受盲件 ID，此时对同件全部四个拟合视图汇报原生残差。主评分用 96²，32² 同口径的材料度量另存。材料评分沿用 SRAP；整体、内部和空气误差同除相同的 \(\lVert m_{\mathrm{ref}}\rVert_2\)，在预测和参考上共用参考半峰阈值计算 IoU。批量输入是带 `predictions` 列表的 JSON，每条提供 `case`、`arm`、`path`，以及可选 `array`、`index`。路径相对列表所在目录解释。迭代求解失败时，在同一保存的 FFT 模型上组装原生有源格 LU 并复核状态缺陷，记录失败迭代和回退费用。仍未解决的复场项目与已算的材料评分同时保留。 [原生正演](physics.py) · [独立真值验证](verify_truth.py) · [评分](scoring.py)

本轮生成用 300 个原生物理右端、12,909 个卷积列，墙钟 6.29 秒。Arb 独立核对全部 **48/48** 份原生误差球，最小严格余量 \(4.10887535\times10^{-6}\)，耗时 94.43 秒。 [无标签的获取和证明统计](ACQUISITION_RECEIPT.json)

## 文件与运行

默认数据根为项目 `.artifacts/z_alg_xudong_cie/b3_holdout`，可用 `--root` 改到同一物理合同的本地路径。`public/` 有两种模型尺寸、公开观测和匿名清单。`offline/` 包含真材、隐藏端口和生成记录。公开加载器只读 `public/`，盲评分需要事先冻结。图示使用事先固定标识的两例，配对汇总按同件、同阵列、同噪声比较，不把四份观测当独立对象。

```sh
python -B research/z_alg_xudong_cie/holdout/data.py --list
python -B research/z_alg_xudong_cie/holdout/generation.py --freeze-record /path/to/freeze.json
python -B research/z_alg_xudong_cie/holdout/verify_truth.py --freeze-record /path/to/freeze.json
python -B research/z_alg_xudong_cie/holdout/scoring.py \
  --freeze-record /path/to/freeze.json --predictions /path/to/prediction_manifest.json \
  --output /path/to/scores.json
python -B research/z_alg_xudong_cie/holdout/summarize.py /path/to/scores.json \
  --baseline dual --output /path/to/paired_summary.json --plot /path/to/paired_differences
python -B research/z_alg_xudong_cie/holdout/materials_figure.py /path/to/scores.json \
  --freeze-record /path/to/freeze.json --output /path/to/fixed_examples
```

运行需 NumPy、SciPy、python-flint、threadpoolctl；出图需 Matplotlib。项目历史脚本与归档物理数组位于原目录。生成、评分和出图以实际冻结的文件为入口。

本轮联合重放可运行 `offline_score.py --wait`。每件只在其全部预定算法输出封存后评分。整个队列完成前，只写本地离线断点；最终生成 `benchmark/synthetic_scores.json` 与 `benchmark/chamber_scores.json`。
