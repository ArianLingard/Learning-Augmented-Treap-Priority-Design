---
author: Arian Lin
---

**函数构造**

根据分桶策略，函数均满足 $priority\mathrm{=}Q\mathrm{(}f(w)\mathrm{)+}U(\mathrm{0},\delta)$，其中 $\delta$ 满足自适应性扰动，因此唯一需要设计的就是 $f(w)$，以满足上述性质和修正策略。

6.1 Softplus平滑双对数函数

为了在消除右端奇点的同时，避免引入如 $ϵ$ 这样缺乏物理意义的经验常数，我们可以构造如下函数：

$$
f_\mathrm{1}(w)\mathrm{=-}\mathrm{log}_\mathrm{2}(\mathrm{log}_\mathrm{2}\frac{\mathrm{1}}{w}\mathrm{+}D(w))
$$

（其中 $D(w)\operatorname{=}{\ln{(1+e^{\gamma w})}-\ln{2}}$，阻尼系数需满足 $\ln(2e-1)\approx1.49<\gamma\lesssim1.71$：下界保证 $D(1)>1$ 以消除右端奇点，上界 $\gamma e^{\gamma}/(1+e^{\gamma})<1/\ln2$ 保证 $f$ 随权重单调不减）

构造逻辑：引入$D(w)=Softplus(\gamma w)-\ln{2}$ 作为平滑正则化阻尼项，实现对数弹性有界。

优势：当 $w\mathrm{\to0}$时项趋于零，函数退化为标准双对数，保留长尾区分度；当 $w\mathrm{\to1}$时提供平滑支撑，消除奇点。

缺陷：仍依赖绝对权重 $w$ 的量纲，不满足尺度不变性。

6.2 线性秩映射函数

$$
f_\mathrm{2}(w)\mathrm{=}K\mathrm{-}⌊\frac{R(w)\mathrm{⋅}K}{n}⌋
$$

(其中 $R(w)\mathrm{\in}[\mathrm{1},n]$为 $w$ 的降序排名)

**构造逻辑**：彻底抛弃权重幅度，只保留降序秩 $R(w)$，线性铺到 $[\mathrm{0},K]$。

**优势**：完全尺度不变；自带离散化；预测器量纲漂移时拓扑不变。

**缺陷**：丢失频率幅度；热门与次热门可能只差一层；退化为近似均匀分层，丢失静态最优性。

6.3 基于线性秩的平滑双对数函数

$$
f_\mathrm{3}(w)\mathrm{=-}\mathrm{log}_\mathrm{2}(\mathrm{log}_\mathrm{2}\frac{\mathrm{1}}{p}\mathrm{+}\epsilon)\text{, }p\mathrm{=}\frac{n-R(w)+1}{n}
$$

其中 $R(w)\mathrm{\in}[\mathrm{1},n]$为 $w$ 的降序排名，$\epsilon$ 为正则常数（如 $1e-9$）。

**构造逻辑**：把秩归一化成伪概率 $p$，套回双对数，恢复长尾层级区分；对数内部加极小常数 $\epsilon$ 消除边界奇点。

**优势**：尺度不变 + 双对数分层；比线性秩更接近 Chen 的拓扑形状。

**缺陷**：本质假设“秩均匀分位”，丢失真实频率密度；静态最优性只保留一部分；精确全局秩需离线或近似维护。

6.4 基于累计分布的平滑双对数函数

**设当前预测器对全量数据的输出分布拟合出一个平滑的累积分布函数（CDF） $F(w)$**

$$
f_\mathrm{4}(w)\mathrm{=-}\mathrm{log}_\mathrm{2}(\mathrm{log}_\mathrm{2}(\frac{\mathrm{1}}{F(w)+1/n})\mathrm{+}\epsilon)
$$

（其中 $\epsilon$ 为正则常数，$+1/n$ 修正保证最轻键 $F=0$ 时对数有定义）

构造逻辑：把离散秩升级为连续 CDF 分位数 $F(w)$，保留预测器输出的相对密度形状；对数内部加极小常数 $\epsilon$ 消除边界奇点。

优势：尺度不变；比秩双对数保留更多局部密度；热点扎堆可被 CDF 斜率自然铺开。

缺陷：$F(w)$ 是估计量，CDF 误差直接传导到层级边界；流式维护成本高；非真实概率，绝对频率幅度丢失。
