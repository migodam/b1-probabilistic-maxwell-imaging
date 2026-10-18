# P1 材料模型、概率测度及局限

## 实际物理对象

采用原24×24 scalar TM/VIE网格，外场与接收几何、波数2pi和已知均匀宿主不变。已知场景支撑包络为r≤.38，包含256个网格单元；它是成像问题的已知ROI，内部材料形状未知，不是由GT或A17成本选出的端口。

宿主chi_b=.1+.04i，epsilon_r=1+chi。设混合分数phi∈[0,1]，材料为

\[
\chi=\chi_b+\phi\,(a+ib),\qquad a\sim U[.8,1.8],\quad b\sim U[.03,.18].
\]

这定义了被动有效标量介质；没有声称任意实际微结构都严格服从同一混合律。phi=0为已知宿主，phi=1为纯包含相；宿主不是空气。材料实部/损耗始终使用共同的复材料，而非每source独立拟合。

## 相关场生成器

z0…z63独立N(0,1)。z0只控制DC偏置，余63个是8×8 DCT频率中去掉DC的空间模式，没有重复DC或共同权重缩放。物理网格上的模式按单位全网格RMS归一化，谱系数与[1+(pi ell)^2(k²+l²)]^{-1}成正比，再统一归一化；ell=.08为有限谱rolloff参数，不是任意形状的硬最小特征宽度，也不声称精确无限域stationary Matérn covariance。

\[
f=.65z_0+\sum_{(k,l)\ne(0,0)} A_{kl}z_{kl},\qquad
\phi=1_{r\le .38}\,h(f/.5).
\]

h(s)=0（s≤−1）、1（s≥1）；中间令t=(s+1)/2，h=6t^5−15t^4+10t^3。它是C2饱和函数，可自然形成纯宿主区、纯材料区和混合过渡，不预先指定孔洞或component数量。解析导数h'=15t²(1−t)²（饱和区为0）。phi数值最大值为1+7e−16，属于浮点舍入，采用1e−12容差；没有物理上放大填充分数。

随机场可以产生低于任务尺度的小岛；相关谱不是硬最小feature约束。任务评价另外要求至少4个连通变化单元及平滑fraction RMS，不能用一个像素声称发现新结构。

## 正确的latent prior与posterior

q=(z64,u2)，u独立U[0,1]，分别线性映射到a和b。因此相对于dz du：

\[
\log p_0(q)=-\tfrac12\|z\|^2-32\log(2\pi),\quad u\in[0,1]^2.
\]

这是正规化的Gaussian×uniform先验，材料law是其decoder推前。饱和产生多对一映射及可能的纯宿主材料原子，不能将其改称材料Lebesgue密度、附加Gram determinant，或把不同latent误当不同材料。没有继承P0归一化Gaussian的weight-ratio law。

局部提案采用z的pCN与u的固定对角反射Gaussian乘积，对原prior可逆；全局对照为独立prior draw。使用完整Maxwell likelihood时，MH比率只剩exp(Phi_old−Phi_new)。这只是核的数学接口，当前128步练习是prior-only，不是posterior采样或混合认证。

## 测量与noise

固定8sources×32complex receivers，forward为总scattered field S j，gain1；逐source打包Re32后Im32。每个complex通道噪声CN(0,sigma²)，实虚方差各sigma²/2。sigma=.00119211798514519，只读取F3已暴露weak观测metadata的sigma标量，没有读取其GT。对P1父对象不重新设20dB；实际SNR需等full-forward后报告。

\[
\Phi(q)=\|F(g(q))-y\|^2/\sigma^2,\quad
\log p(y|q)=-256\log(\pi\sigma^2)-\Phi(q).
\]

posterior相对于同一latent prior定义，不profile gain，不把noise-ball或DeltaPhi≤2当硬支持。实际y尚未生成，evidence/event posterior mass/CDF均NOT_RUN。

## 事件与材料度量

全材料RMS=||delta chi||/24，主任务s=.02、次任务.04，均在样本输出和任何Maxwell结果前冻结。结构检查为fraction改变≥.15的4连通块至少4单元，及1pixel平滑fraction RMS≥.01。

事件只在输出后由phi≥.25/.5/.75计数：foreground4连通、host8连通、至少4单元。五个主要事件为无足够大threshold component、多component、单component无闭合宿主区、一个闭合宿主区、多个闭合宿主区。这是阈值材料事件，不是generator输入；不能把threshold pocket等同严格空气孔。

当前模型限定已知ROI、均匀host和二维TM；没有一般未知三维几何、色散、仪器nuisance或真实实验设备验证。更丰富先验并不保证观测有多个重要模式。
