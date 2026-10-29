# 数字8：新非参数全波探索

只对数字8生成了真实PDE观测并进行了实际110实变量反演；1..9图片是公共字形定义，不能称为九个成像结果。该探索与原椭圆盲测分开。

未知量为墙6个实部/6个loss，以及7×7目标全部49个contrast和49个loss。无模板分类、无ellipse、未提供笔画support；目标为独立box且无空间正则。笔画宽0.2λ₀，对象1.0×1.4λ₀，两个真孔约0.6×0.4λ₀。模型材料格0.2λ₀；FEM公共targetstep0.1λ₀、P3，生成与反演同离散，不是独立细网格PASS。

完整110变量weighted-adjoint L-BFGS-B，每次完整梯度3FE/24RHS，至多两start×40iter，真实calls/450FE/3600RHS/18:30先到停止。初始化为全target统一contrast0.1/loss0.03，不知道数字ID或真墙。原频率、四发十六收和每实虚坐标σ0.01不变。

结果图truth_vs_reconstruction.png的contrast色轴固定0..8、loss固定0..0.25；不得自动拉色轴制造成功。support示例χ≥1；以下同时报告0.5/1/2阈值，不选择最好者，也不新增PASS门。

```json
{
  "stops": [
    {
      "start": 0,
      "calls": 44,
      "stop": "STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT"
    },
    {
      "start": 1,
      "calls": 44,
      "stop": "STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT"
    }
  ],
  "metrics": {
    "target_complex_contrast_relative_L2": 0.779265072370431,
    "target_epsilon_complex_relative_L2_f1": 0.5791412892823169,
    "contrast_RMSE": 1.738434759077241,
    "loss_RMSE": 0.04528445362387194,
    "contrast_centroid_error": 0.04803406002566643,
    "offstroke_contrast_mass_fraction": 0.42267348361767854,
    "padding_contrast_mass_fraction": 0.2672762609367155,
    "truehole_contrast_mass_fraction": 0.15539722268096298,
    "wall_epsilon_RMSE": 0.09203300814272711,
    "wall_loss_RMSE": 0.07063163448153766,
    "observation_white_norm": 19.205094524468976,
    "observation_white_RMS": 0.9800558759764195,
    "same_discrete_clean_prediction_white_RMS": 0.4559952093726813,
    "threshold_diagnostics": [
      {
        "threshold": 0.5,
        "support_IoU": 0.47368421052631576,
        "estimated_enclosed_holes": 4,
        "true_hole_representative_components": [
          0,
          6
        ],
        "two_true_holes_remain_distinct_enclosed_background": false
      },
      {
        "threshold": 1.0,
        "support_IoU": 0.40540540540540543,
        "estimated_enclosed_holes": 1,
        "true_hole_representative_components": [
          0,
          4
        ],
        "two_true_holes_remain_distinct_enclosed_background": false
      },
      {
        "threshold": 2.0,
        "support_IoU": 0.35714285714285715,
        "estimated_enclosed_holes": 0,
        "true_hole_representative_components": [
          2,
          2
        ],
        "two_true_holes_remain_distinct_enclosed_background": false
      }
    ]
  },
  "cost": {
    "full_fe": 276,
    "rhs": 2172
  }
}
```

以上属于限定时间、粗材料格、单个预声明digit8的实际探索。未收敛/失败不得解释为物理无解；不能证明1..9复杂目标已恢复，也不改变原R3结论。

## 目视结论与执行停止

连续contrast图并未清楚恢复数字8：笔画断裂、空白出现显著材料、两个孔洞未可靠分开。固定阈值0.5/1/2均未保留两个真孔中心对应的独立闭背景。不能将本次实际反演输出称为成功识别8。两个start均达到40iter上限，故也不能据此宣布该观测物理上无法成像数字。

唯一直方mixed方向梯度相对差为2.64e-8，只核对实现导数，不是成像证据。实际276FE/2172RHS，online124.04CPU秒/127.32wall秒；生成2.92CPU秒。18:15:09UTC已封存，随后仅评分/文件工作，无新增其他数字PDE。见RUN_DIGITS.md。

## 0PDE：真值与封存候选的同一训练目标对照

真值恰在相同7×7合法材料类中，使用已存clean均值、y和相同墙prior即可精确计算该离散真值目标，无需新PDE。各start点只按既有训练目标最小规则从已封存候选读出，不根据真值选点。

```json
{
  "time_unix": 1793298158.011735,
  "new_PDE": 0,
  "true_model_in_exact_same_legal_7x7_class": true,
  "truth_full_objective": 213.69369246714248,
  "truth_data_objective": 212.60283758158374,
  "truth_wall_prior_objective": 1.090854885558739,
  "actual_noise_raw_complex_norm": 0.2062051588014149,
  "actual_noise_white_norm": 20.62051588014149,
  "n_complex_channels": 192,
  "sigma_each_real_imag": 0.01,
  "per_start_saved_points": [
    {
      "start": 0,
      "selected_saved_call": 44,
      "selection": "minimum already-sealed full measured objective within start; no truth selection/no refit",
      "trace_objective": 185.70888701976637,
      "recomputed_full_objective": 185.70888701976637,
      "data_objective": 184.4178278468941,
      "wall_prior_objective": 1.2910591728722627,
      "full_objective_minus_truth": -27.984805447376118,
      "lower_training_objective_than_truth": true,
      "clean_mean_distance_raw_complex_norm": 0.08935644704933203,
      "clean_mean_distance_white_norm": 8.935644704933203,
      "distance_not_divided_by_channel_count": true,
      "last_saved_call": 44
    },
    {
      "start": 1,
      "selected_saved_call": 44,
      "selection": "minimum already-sealed full measured objective within start; no truth selection/no refit",
      "trace_objective": 218.54041899452056,
      "recomputed_full_objective": 218.54041899452056,
      "data_objective": 217.33642135250597,
      "wall_prior_objective": 1.2039976420145988,
      "full_objective_minus_truth": 4.846726527378081,
      "lower_training_objective_than_truth": false,
      "clean_mean_distance_raw_complex_norm": 0.1005987049247578,
      "clean_mean_distance_white_norm": 10.05987049247578,
      "distance_not_divided_by_channel_count": true,
      "last_saved_call": 44
    }
  ],
  "limits": "Finite noisy training-objective comparison only; not exact fiber, global nonidentifiability, global optimum, or physical impossibility proof; no online change or hyperparameter adjustment"
}
```

错误材料图即使训练目标低于真值，也不能据较低residual宣称形状恢复。这里的clean均值距离是完整192复通道Euclidean/white范数，不除通道数。此有限带噪对照不证明exact fiber、全局不可辨或物理无解，也未触发调参/重选/重跑。
