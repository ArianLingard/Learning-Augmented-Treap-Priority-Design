在 Treap 的大根堆中，priority 越大越靠近根。Priority
函数的目的就是实现"高频访问的节点更浅"。

# 宏观骨架构造

有一权重集合是
$W_{n} = \left\{ w_{1},w_{2},\ldots,w_{n} \right\}$，存在一个映射函数
$f:\mathcal{D}\mathbb{\rightarrow R}$（其中 $\mathcal{D}$
为权重的取值域，记
$d_{\min} = \inf\mathcal{D}$，$d_{\max} = \sup\mathcal{D}$），满足以下性质：

## 性质 1：严格单调保序性

数学表述：设 $w_{x} > w_{y}$，则
$\mathbb{E}\left\lbrack P(x) \right\rbrack\mathbb{> E}\left\lbrack P(y) \right\rbrack\text{\:\,} \Longleftrightarrow \text{\:\,}w_{x} > w_{y}$

工程目的：系统成立的基础，确保热数据比冷数据更靠近根。

## 性质 2：左端点导数发散

数学表述：在 $w \rightarrow d_{\min}^{+}$时，
$\mid f'(w) \mid \rightarrow \infty$

工程目的：放大长尾极小区间内的微小差异，防止海量冷门节点在浮点数中坍缩成同一个值。

## 性质 3：右端点导数有界

数学表述：在
$w \rightarrow d_{\max}^{-}$时，$\underset{w \rightarrow d_{\max}^{-}}{lim\, sup}f'(w) \leq L_{H} < \infty$

工程目的：抑制超级热点，防止绝对热点独占根节点导致顶部局部退化。

## 性质 4：全局极差约束

数学表述：$M(n) = \max_{w \in W_{n}}f(w) - \min_{w \in W_{n}}f(w) = \Theta\left( \log n \right)$

工程目的：约束极差的大小，确保查询的最坏比较次数被压缩在
$O\left( \log n \right)$量级；同时防止浮点数出现上溢或下溢，维持数值表示的有效性。

# 局部修正补丁

## 修正 1：针对冷点差值过小精度截断的分桶策略

性质2貌似已经完备，但是在实际工程落地时存在缺陷，理由如下：

数学上，两点间的差值可以近似为：$\Delta f \approx f'(w) \cdot \Delta w$

当 $w \rightarrow 0$时，性质 2 保证了
$f'(w) \rightarrow \infty$。但是，如果 $\Delta w$ 趋于 0 的速度比
$f'(w)$趋于无穷大的速度更快，那么它们的乘积 $\Delta f$
依然可以是一个极小的数。

计算机是 IEEE 754
浮点数系统，精度是有限的，在数学上可以被区分的两点，在计算机内也许无法被区分。

因此，我们提出了分桶策略，具体操作如下：

### 离散化

数学表述：取离散化算子 $Q( \cdot )$，使得 $priority$ 函数具有
$Q\left( f(w) \right)$ 的结构。

工程目的：强制剥离计算机有限精度下的微小预测差异，将无限细分的实数域量化。确保同层内节点在
IEEE 754 浮点数表示下具有至少 $\Omega(1)$ 的绝对差值，防止因精度截断导致
priority 相等，从而避免堆性质失效和查询复杂度退化为 $O(n)$。

### 微观扰动

为了使得同层结点之间亦有区分度，引入 $\varepsilon$ 作为微观扰动。

在 Treap 中，节点 $x$
的深度等于它的祖先个数。所以，期望深度就是所有其他节点 $y$ 成为 $x$
祖先的概率之和：$\mathbb{E}\left\lbrack \text{depth}(x) \right\rbrack = \sum_{y \neq x}^{}P\left( y\text{是}x\text{的祖先} \right)$。

根据 Treap 的堆性质：$y$ 是 $x$ 的祖先，当且仅当 $y$ 的 priority 是区间
$\left\lbrack x,y \right\rbrack$ 内所有节点中最大的。

为了使得查找期望为 $O\left( \log n \right)$，选择均匀随机分布，即取
$\varepsilon = U(0,\delta)$。

此时期望为：$\mathbb{E}\left\lbrack \text{depth}(x) \right\rbrack = \sum_{y \neq x}^{}\frac{1}{\mid x - y + 1 \mid} = \sum_{k = 1}^{n}\frac{1}{k} = O\left( \log n \right)$。

### 跃迁隔离

在分桶之后，性质 1 即被弱化成了"跨层严格保序，层内随机"。

为了防止高频节点与低频节点在堆结构中的相对位置因预测噪声而发生逆转，必须实现"绝对的层级隔离"。

设高频节点集合为 $H$，低频节点集合为 $L$，且对任意 $x \in H,y \in L$，有
$w_{x} > w_{y}$。

映射后的期望隔离下界为：

$$\inf_{x \in H}\mathbb{E}\left\lbrack f\left( w_{x} \right) \right\rbrack - \sup_{y \in L}\mathbb{E}\left\lbrack f\left( w_{y} \right) \right\rbrack \geq \Delta > 0$$

在真实系统里，由于预测器的噪声或浮点数的截断误差，基础分会发生抖动若上下两层过于贴近，导致精度截断。

为了确保跨层不发生翻转，例如，可以取
$\delta = \Delta/2$、$\delta = \sqrt{\Delta}$、$\delta = log(\Delta)$
等。

### 最高频结点边界修正

由于注入了随机扰动
$\epsilon \sim U\left( 0,\delta \right)$，因此还要求最高频节点的基础分加上最大扰动后，依然不能超过
$O\left( \log n \right)$的量级（设高频节点集合为 $H$），即：

$$\sup_{w \in W_{H}}f(w) + \delta \leq O\left( \log n \right)$$

### 针对性质1的修正

**跨层严格保序**：若
$Q\left( f\left( w_{x} \right) \right) > Q\left( f\left( w_{y} \right) \right) + \delta_{\max}$，则
$P(x) > P(y)$ 必定成立；

**层内随机**：同层只保证期望/概率意义下的优先级，不保证任意两点的严格序。

### 针对性质4的修正

设离散化后的层号为 $\tau(w) = Q\left( f(w) \right)$，层间距为
$L$，则要求：

$$\max_{w \in W_{n}}\tau(w) - \min_{w \in W_{n}}\tau(w) = O\left( {\log\log}n \right)$$

因而 priority 的数值极差满足

$$R(n) = {\max{priority}}(w) - {\min{priority}}(w) \leq LT(n) + \gamma = O\left( {\log\log}n \right).$$

### 分桶策略总结

$$priority = Q(f(w)) + U\left( 0,\delta \right)$$

分桶策略通过强加一个离散化操作，把差异微小的数字都划分在 $w$
发生极微小变化时，强行让 $f(w)$ 发生至少 $\Omega(1)$ 的阶跃。

这样一来，可以确保不同层的节点在 IEEE 754 浮点数表示下具有至少
$\Omega(1)$ 的绝对差值，防止因精度截断导致 priority
相等，从而避免堆性质失效和查询复杂度退化为 $O(n)$。

分桶策略把真实世界的$\Delta w$，转变为计算机里的$\Delta w$。因为真实数据之间的差距计算机也许无法识别，但是随机数产出的结果计算机一定是可以识别的。

## 修正2：针对热点潜在相对误差过大的修正策略

性质3貌似已经完备，但是在实际工程落地时存在缺陷，理由如下：

数学上，两点间的差值可以近似为：$\Delta f \approx f'(w) \cdot \Delta w$

当 $w \rightarrow d_{\max}^{-}$ 时，性质 3 保证了 $f'(w)$ 为一常数。但
$f'(w)$ 过大时，仍然可能导致乘积 $\Delta f$ 依然可以是一个极大的数。

在实际工程中，若 $priority$ 函数设计得不好，一个极小的 KL
散度扰动，会导致节点的 priority 发生剧烈跳跃。（KL 散度衡量的是预测分布
$Q$和真实分布 $P$之间的信息差异。）

> *"We assume wlog that* $\tau_{x} \geq 0$*for any* $x$*. Otherwise,*
> $\tau_{x} < 0$*implies* $w_{x} = \Omega(1)$*, which can hold for only
> a constant number of items. So, we can always put them at the top of
> the Treap, which increases the depths of other items by a constant."\
> ------《On the Power of Learning-Augmented Search Trees》Theorem 2.4*
>
> （我们不失一般性地假设 Tier $\tau_{x} \geq 0$。如果
> $\tau_{x} < 0$，意味着权重
> $w_{x}$非常大，但这最多只有常数个元素。所以我们直接把它们放到 Treap
> 的最顶层，这只会让其他元素的深度增加一个常数。）

该论文中的函数（$f(w) = - \log_{2}\log_{2}\left( \frac{1}{w} \right)$）也存在这一问题，但作者[在理论推导时，直接认定超级热点最多就几个，直接把它们放在最上面当根节点，忽略它们导致的数值溢出。]{.smallcaps}

现在笔者给出其他的应对策略：

### 平滑正则化

在函数的内部（尤其是分母或对数内部）直接加一个极小常数
$\epsilon$。例如：

$$f_{\text{smooth}}(w) = - \log_{2}\log_{2}\left( \frac{1}{w + \epsilon} \right)$$

取 $w' = w + \epsilon$，强制让权重永远不真正等于 0 或 1。

数学本质：通过改变函数的定义域或修改内部结构，消除函数在边界处的奇点（Singularity），使得函数变得有界，导数变得有限。

代价是引入了无法用 KL
散度理论解释的系统偏差，并破坏了函数在边界处的光滑性与单调性，只是修正函数右侧的补丁。

### 对数弹性有界

在现实系统中，预测器输出的权重
$w_{x}$（比如用户点击某个词的概率）往往不是完美的。它最常犯的错误不是"加减法"的错误，而是乘性误差。

一旦发生这种乘性缩放，节点之间的相对距离就会彻底畸变。原本距离很近的两个节点，可能瞬间跨越好几个层级；原本距离很远的节点，可能被压缩在一起。

为了严格刻画"比例缩放"，经济学中有一个概念是点弹性，定义为：

$$E(w) = \frac{\Delta f/f}{\Delta w/w} = \frac{df}{dw} \cdot \frac{w}{f(w)} = \frac{d\ln f(w)}{d\ln w}$$

弹性描述了当输入 $w$变化 $1\%$时，输出 $f(w)$变化的百分比。

如果我们要求对数弹性有界（$\mid E(w) \mid \leq E$），意味着函数
$f(w)$不能发生"指数级爆炸"，它最多只能是幂函数级别的增长（$f(w) \sim w^{\alpha}$，其中
$\alpha \leq E$）

这样一来，无论预测器把 $w$ 放大了多少倍，$f(w)$ 的增长率都不会产生畸变。

数学表述：存在常数 $E > 0$，使得在
$w \rightarrow d_{\max}^{-}$区域，满足：
$\left. \mid\frac{d\ln\left| f(w) \right|}{d\ln w} \right.\mid \leq E$

### 概率校准

现代机器学习模型（特别是基于 Softmax 的神经网络或大型
RNN）输出的概率，并不是真实的概率，而是置信度得分。

概率校准的目标，是寻找一个严格单调递增的映射函数
$C( \cdot )$，将模型的原始输出 $q$，转化为真实的经验概率 $\widehat{p}$：

$$\widehat{p} = C(q)$$

使得对于所有的 $q$，预测为 $q$的样本中，真实正样本的比例恰好是 $C(q)$。

常用的工程手段是温度缩放（Temperature Scaling）：

$${\widehat{p}}_{i} = \frac{\exp\left( z_{i}/T \right)}{\sum_{j}^{}\exp\left( z_{j}/T \right)}$$

其中 $z_{i}$是模型 logits 输出，$T$ 是温度参数（$T > 1$
时，分布被平滑；$T < 1$
时，分布被尖锐化）。通常在验证集上通过最小化负对数似然（NLL）来求得最优的
$T$。

在进入映射函数前，先用温度缩放（Temperature
Scaling）把模型过度自信的输出平滑回合理区间。具体的工程逻辑与priority函数设计无关，放在预测模型的设计中详细深入。

## 修正3：针对预测器模型变化的策略

### 尺度序列不变性

**解决预测器模型发生变化的问题（包括权重值域的变化 +
对权重排序方式的变化）**

问题一：预测模型输出的权重可能是概率（0\~1）、相对频率（1\~n）或者某种得分，同一批数据算出来的
Tier 会完全不同。Treap的结构会发生变化。

问题二：如果模型对权重做了任意严格单调递增的变换 $g( \cdot )$（例如
$w^{2}$或者 $e^{w}$），Treap的结构会发生变化。

这意味着 $f(w)$和
$f\left( g(w) \right)$必须在拓扑上保持一致，通过设计一种基于局部秩或分位数的保序映射，可以绕过对全局归一化的依赖。

定义秩函数 $R:W \rightarrow \left\{ 0,1,\ldots,n - 1 \right\}$。映射函数
$f:W\mathbb{\rightarrow R}$满足：对于任意严格单调递增变换
$g( \cdot )$，有 $R(w) = R\left( g(w) \right)$，从而 $f(w)$**与**
$f\left( g(w) \right)$**的秩等价**。

工程目的：确保priority的映射结果仅依赖于权重的相对排序，而不受全局量纲缩放的影响。允许在不重建树结构的前提下，动态替换预测器或进行数据流的全局归一化。

# Chen函数的缺陷

$$priority = - \left\lfloor \log_{2}\log_{2}\left( \frac{1}{w} \right) \right\rfloor\  + U\left( 0,1 \right)$$

## 无法真正解决"热点潜在相对误差过大"的问题

**数学表现**：当
$w \rightarrow 1^{-}$时，$\log_{2}(1/w) \rightarrow 0^{+}$，外层对数导致
$f(w) \rightarrow + \infty$，其一阶导数 $f'(w) \rightarrow \infty$。

**工程后果**：无法从解析式上限制超级热点的优先级绝对峰值。当预测器给出极为自信的输出（如
$w = 0.9999$）时，映射值会趋近于无界，极易引发浮点数溢出（Infinity）。

**原论文处理**：原论文未在函数层面进行修正，而是通过引入强假设（WLOG）予以规避：在
Theorem 2.4
证明中假设"超级热点仅存在常数个"，直接将其置于树顶，从而忽略该区间的数值发散问题。

## 工程层面极差约束失效

**数学表现**：在理想状态下，双对数映射可将值域极差压缩在
$O\left( {\log\log}n \right)$内。但由于右侧奇点的存在，一旦输入逼近
1，函数值失控。极差 $M$ 根本没有绑定数据规模
$n$，而是被预测器的绝对数值直接控制。如果 $w_{x}$ 没有被归一化为
$[0,1]$，极差会很大。导致树的期望高度过高。

**工程后果**：理论极差界在真实机器学习输出下失效，无法保证防溢出的物理底线。

## 不满足尺度序列不变性

**数学表现**：映射函数强依赖于绝对数值量纲 $w$。对于任意严格单调递增变换
$g(w) = c \cdot w$，由于双对数的非线性平移效应，$f(w) \neq f\left( g(w) \right) + C$。

假设预测器给出的权重发生缩放
$g(w) = c \cdot w$。$\log\log\left( 1/\left( c \cdot w_{x} \right) \right) = \log\left( \log\left( 1/w_{x} \right) - \log c \right)$。因为外层是双对数，这个平移项
$- \log c$无法被轻易提取出来，它会导致**曲线发生严重的形变**。如果换个数据集，预测器改变，所有节点的相对层级分布会彻底重排，原本在顶层的可能掉到下面。它根本不具备真正的"尺度不变性"。

**工程后果**：当预测器模型发生更替（如从输出概率变更为未归一化频数）或数据流发生量纲漂移时，所有节点的
Tier
分布会发生全局畸变与偏移。系统无法做到仅依赖"相对秩"来维持拓扑结构，必须进行全量重建。

## 不满足微观扰动隔离

**表现**：原论文设定层间距 $\Delta = 1$，并引入均匀扰动
$\delta \sim U\left( 0,1 \right)$。

**工程后果**：扰动幅度的上界恰好等于层间距，未预留安全缓冲区。在浮点数截断误差或预测噪声的干扰下，极易发生跨层边界碰撞（即高频层的最低分与低频层的最高分发生翻转），缺乏绝对的层级隔离鲁棒性。

# 构造函数

根据分桶策略，函数均满足
$priority = Q(f(w)) + U\left( 0,\delta \right)$，因此唯一需要设计的就是
$f(w)$，以满足上述性质和修正策略。

## 平滑双对数函数

$$f_{1}(w) = - \log_{2}\left( \log_{2}\frac{1}{w} + \epsilon \right)$$

**构造逻辑**：保留 Chen 双对数骨架，只在内部加 $\epsilon$把
$w \rightarrow 1$的奇点抬离原点。

**优势**：改动最小；左端长尾分层保留；右端有界，消除浮点上溢。

**缺陷**：仍依赖绝对权重 $w$ 的量纲；$\epsilon$
无理论取值；不满足尺度不变性。

## 线性秩映射函数

$$f_{2}(w) = K - \left\lfloor \frac{R(w) \cdot K}{n} \right\rfloor$$
(其中 $R(w) \in \left\lbrack 1,n \right\rbrack$为 $w$ 的降序排名)

**构造逻辑**：彻底抛弃权重幅度，只保留降序秩 $R(w)$，线性铺到
$\left\lbrack 0,K \right\rbrack$。

**优势**：完全尺度不变；自带离散化；预测器量纲漂移时拓扑不变。

**缺陷**：丢失频率幅度；热门与次热门可能只差一层；退化为近似均匀分层，丢失静态最优性。

## 秩双对数映射函数

$$f_{3}(w) = - \log_{2}\left( \log_{2}\frac{1}{p} + \epsilon \right),\ p = \frac{n + 1 - R(w)}{n + 1}$$
(其中 $R(w) \in \left\lbrack 1,n \right\rbrack$为 $w$ 的降序排名)

**构造逻辑**：仍只信秩，但把秩归一化成伪概率
$p$，再套回双对数，恢复长尾层级区分。

**优势**：尺度不变 + 双对数分层；比线性秩更接近 Chen 的拓扑形状。

**缺陷**：本质假设"秩均匀分位"，丢失真实频率密度；静态最优性只保留一部分；精确全局秩需离线或近似维护。

## 基于概率积分的"平滑秩"双对数

设当前预测器对全量数据的输出分布拟合出一个平滑的累积分布函数（CDF）$F(w)$。映射函数定义为：

$$f_{4}(w) = - \log\left( \log\left( \frac{1}{F(w) + \frac{1}{n}} \right) + \epsilon \right)$$

**构造逻辑**：把离散秩升级为连续 CDF 分位数
$F(w)$，保留预测器输出的相对密度形状，再用 $+ 1/n$避开右端奇点。

**优势**：尺度不变；比秩双对数保留更多局部密度；热点扎堆可被 CDF
斜率自然铺开。

**缺陷**：$F(w)$ 是估计量，CDF
误差直接传导到层级边界；流式维护成本高；仍非真实概率，绝对频率幅度丢失。

# 构造函数的缺陷

## 绝对频率静态最优

给定访问序列 $X = \left( x,(1),\ldots x(m) \right)$，令 $f_{x}$为键
$x$的出现次数，$p_{x} = f_{x}/m$。若数据结构处理 $X$的总访问代价为

$$O\text{ }\left( \sum_{x}^{}f_{x}\log_{2}\frac{m}{f_{x}} \right) = O\text{ }\left( m \cdot Ent(p) \right)$$

则称其在常数因子内达到**静态最优**。该界是静态最优 BST 的熵下界，Chen et
al. 通过设 $w_{x} = p_{x}$达到。

## 方法对比

Chen 的复合 priority 在 $w_{x} = p_{x}$ 时给出期望深度：

$$O\left( \log_{2}\left( \frac{1}{p_{x}} \right) \right)$$

因而总代价 $O\left( m \cdot Ent(p) \right)$
在常数因子内达到绝对频率静态最优。

本文的
$f_{1}$保留绝对权重，因此在预测准确时继承该界；$f_{2},f_{3},f_{4}$
仅依赖秩或分位数，达到的是**分位数静态最优**：存在由秩/CDF 诱导的伪分布
$\widetilde{p}$，使总代价为

$$O\left( \sum_{x} f_x \log_2 \frac{1}{\tilde{p}_x} \right) = O\left( m \cdot Ent(p \parallel \tilde{p}) \right)$$

它与绝对频率静态最优的差为：

$$O\text{ }\left( m \cdot D_{KL}\left( p \parallel \widetilde{p} \right) \right) = O\text{ }\left( \sum_{x}^{}f_{x}\log_{2}\frac{p_{x}}{{\widetilde{p}}_{x}} \right)$$

本文构造的函数达到伪分布下的静态最优，损失由
$D_{KL}\left( p \parallel \widetilde{p} \right)$控制。

### 线性秩

设降序秩
$R_{x} \in \left\lbrack 1,n \right\rbrack$，均匀秩假设下的伪概率

$${\widetilde{p}}_{x} = \frac{n + 1 - R_{x}}{\sum_{y}^{}(\left. \ n + 1 - R_{y} \right)} = \frac{2\left( n + 1 - R_{x} \right)}{n(n + 1)}$$

损失：

$$\sum_{x}^{}f_{x}\log_{2}\frac{p_{x}}{{\widetilde{p}}_{x}} = m \cdot D_{KL}\left( p \parallel \widetilde{p} \right)$$

当真实频率与秩近似线性（如均匀分布）时损失小；当分布极偏斜时损失大。这就是
$f_{3}$在对抗分布下明显变差的原因。

### 秩双对数

同样取
${\widetilde{p}}_{x} = \frac{2\left( n + 1 - R_{x} \right)}{n(n + 1)}$，但代入双对数。深度为
$O\left( \log\left( 1/{\widetilde{p}}_{x} \right) \right)$，损失同上。

### CDF 平滑秩

取
${\widetilde{p}}_{x} = F\left( w_{x} \right) - F\left( w_{x}^{-} \right)$或近似
$F\left( w_{x} \right)$（视定义）。深度
$O\left( \log\left( 1/F\left( w_{x} \right) \right) \right)$。若
$F$是真实频率 CDF，则 $F\left( w_{x} \right)$是尾部概率，一般不等于
$p_{x}$；损失可用这个参数界定：

$$\sum_{x}^{}f_{x}\left. \mid\log_{2}\frac{p_{x}}{F\left( w_{x} \right)} \right.\mid$$

只有当分布重尾、$F\left( w_{x} \right) \approx p_{x}$
时才接近绝对静态最优。

## 测定代价

当预测器输出经过概率校准且保留绝对幅度时，取
$w_{x} = {\widehat{p}}_{x}$的 $f_{1}$在
${\widehat{p}}_{x} \approx p_{x}$下恢复 Chen 的静态最优界。

当预测器只能给出排序或量纲不稳定时，$f_{3},f_{4}$ 以
$m \cdot D_{KL}\left( p \parallel \widetilde{p} \right)$的代价换取尺度不变性；该代价在秩--频率关系接近线性或
CDF
校准良好时较小，在极端偏斜且秩均匀假设失效时增大。实验将报告静态最优比以量化该代价：

$$\rho = \frac{\text{实测总代价}}{OPT}$$

其中 $OPT = \min\sum_{x}^{}f_{x} \cdot depth(x)$ 为最优静态 BST 代价（用 Knuth 动态规划按实际访问频率求得，$depth$ 为节点访问次数、根 = 1），因此 $\rho \geq 1$ 恒成立，越接近 1 越接近静态最优。香农熵 $\sum_{x}^{}f_{x}\log_{2}\left( m/f_{x} \right)$ 是叶子二叉码（按边数、每比较 1 bit）的下界，与内节点树三路比较的口径不同，用它作分母会得到 $\rho < 1$ 的假象。
