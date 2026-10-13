# 连续参考审计

|案例|65/97 CDF最大差|97/自适应CDF最大差|5×5/7×7 response的CDF差|独立事件质量最大差|
|---|---:|---:|---:|---:|
|distinguishable|2.875e-14|1.189e-12|1.992e-04|3.402e-21|
|competitive|2.442e-15|1.776e-15|3.350e-05|4.441e-16|
|information_limited|2.665e-15|8.882e-16|1.001e-06|0.000e+00|
|gain_mismatch|6.706e-14|3.878e-12|2.407e-04|6.573e-14|
|gain_marginalized|3.886e-14|3.009e-14|1.402e-04|3.553e-15|
|representation_mismatch|4.673e-13|4.050e-13|2.546e-04|7.022e-15|

使用129点主CDF、33点独立CDF及直接quantile反解；全部网格差为经验差。主积分25/41/65/97，独立二维分区5/9。5×5/7×7 response差单独记录。四项历史原生full-grid失败保留；四个历史treams原始oracle数组缺失保留。

详细逐事件/坐标、moments、SD、quantile、evidence、CDF误差及adaptive收敛状态见results/SUMMARY.json和cases/*_curves.npz。
