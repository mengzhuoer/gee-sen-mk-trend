# 01 · Sen's slope：趋势幅度

## 1. 它回答什么问题

“这个像元在 2001—2022 年间，**平均每年变化多少**？”

Sen's slope（Theil–Sen estimator，简称 Sen 斜率）不拟合一条让残差平方和最小的直线，
而是把所有年份两两配对，算出每一对之间的斜率，再取这些斜率的**中位数**：

$$
\beta = \mathrm{median}\left\{\frac{x_j-x_i}{j-i}\;|\;1 \le i < j \le n\right\}
$$

其中 $x_i$ 是第 $i$ 年的值，$n$ 是年数。年份不等距时，把 $j-i$ 换成真实时间差即可。

## 2. 为什么用中位数而不是最小二乘

最小二乘（OLS）对异常值敏感：一个极端年份（火灾、干旱、云污染残留）就能把斜率拉偏。
Sen 斜率取中位数，只要少于约 29.3% 的配对斜率被异常值污染，结果就不受影响，
因此适合噪声大、缺测多、年际波动强的遥感序列。

它也不要求残差服从正态分布，这是它比 OLS 更适合遥感长时序的主要原因。
但要注意：**Sen 斜率只估计幅度，不判断显著性**。显著性由 Mann–Kendall 检验给出。

## 3. 配对数量与计算量

| 年数 $n$ | 配对数量 $n(n-1)/2$ |
|---:|---:|
| 10 | 45 |
| 22 | 231 |
| 24 | 276 |

对 $n=22$，每个像元要算 231 个斜率再取中位数。100 万像元就是 2.31 亿个斜率，
所以本地实现要分块（本仓库默认每块 6 万像元），GEE 则用 `ee.Reducer.sensSlope()` 在服务端聚合。

## 4. 结果怎么读

- $\beta > 0$：增加；$\beta < 0$：下降。
- $|\beta|$ 越大，单位时间变化越快。
- 单位 = 变量的单位 / 年：
  - NDVI（无量纲）→ NDVI·a⁻¹；
  - GPP（g C/m²/year）→ g C/m²/yr²。
- 斜率的中位数、分位数（p2/p50/p98）适合在论文里描述区域整体幅度；
  逐像元斜率用于制图。

## 5. GEE 中的写法

GEE 的 `ee.Reducer.sensSlope()` 接收一个双波段影像集合，两个波段分别是
**时间**和**变量值**，输出 `slope` 与 `offset` 两个波段：

```javascript
var withYear = annualCollection.map(function (img) {
  var year = ee.Number(img.get('year'));
  var yearBand = img.select('value').multiply(0).add(year).rename('year');
  return yearBand.addBands(img.select('value'));
});

var sen = withYear
  .select(['year', 'value'])
  .reduce(ee.Reducer.sensSlope());

var slope = sen.select('slope').rename('SenSlope');
```

两个容易踩的点：

1. **年份波段要从变量波段派生**（`value.multiply(0).add(year)`），
   这样年份和变量共享同一网格与掩膜；直接用 `ee.Image.constant(year)` 可能落到默认网格。
2. **`reduce` 前集合要按年份排序**，并保证每个像元年份序列完整；
   本仓库用“完整序列掩膜”（22 年全有效）保证这一点。

## 6. 本地实现

`code/python/sen_mk_core.py`：

```python
i, j = np.triu_indices(n, 1)
dt = (years[j] - years[i])[:, None]
slopes = (x[j] - x[i]) / dt
slope = np.median(slopes, axis=0)
```

`tests/test_algorithms.py` 会把向量化结果与“显式双重循环 + `np.median`”逐像元对比，
确保没有把年份对或时间差写错。

## 7. 局限

- 只描述**单调趋势**；先升后降的序列会被压成一个接近 0 的平均斜率。
- 不能处理周期性（季节已在年度合成中消除，但年际振荡仍在）。
- 斜率大小受起始年/结束年影响：如果首尾恰好是极端年份，幅度会偏大或偏小；
  这是所有趋势估计的共同问题，报告时应给出完整曲线（见案例的年际图）。
- 斜率不为 0 不等于趋势显著；必须结合 MK 的 p 值。
