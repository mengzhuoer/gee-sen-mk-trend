# 04 · 年度指标：Sen + MK 之前最关键的一步

Sen + MK 处理的是“年度序列”。把原始遥感产品合成为一年一个值，这个**年度指标**的选择
往往比检验本身更影响结论。原教程用 NDVI 的“生长季最大值”，把它直接套到 GPP 上是错的。

## 1. 先判断变量是“状态量”还是“通量”

| 类型 | 例子 | 年度指标的物理含义 | 推荐合成 |
|---|---|---|---|
| 状态量 | NDVI、EVI、LAI、地表温度 | 某一时刻/时段的水平 | 生长季**最大值**或**均值** |
| 通量 | GPP、NPP、ET、降水 | 一段时间内的累计量 | 全年/生长季**累加** |
| 频率/事件 | 火灾次数、干旱月数 | 一年内发生次数 | **计数** |

GPP 是光合作用在时间上的累计通量。取最大值等于只统计“光合最强的那 8 天”，
会丢掉季节长度、复种、生长季起止时间等信息；只有全年累加才与“年 GPP”对应。

## 2. NDVI 原版口径（文章）

- 产品：`MODIS/061/MOD13A1`，500 m，16 天，NDVI 比例因子 0.0001；
- 时段：每年 4 月 1 日—11 月 1 日（生长季）；
- QA：保留 `SummaryQA ≤ 1`（0 = good，1 = marginal），剔除冰雪/云；
- 数值范围：NDVI 原始值在 `[-2000, 10000]` 内；
- 年度值：生长季内所有有效观测的 **最大值**；
- 存储：先 `toInt()` 保留原始整数（保证 ties 可精确判断），另建一个
  `× 0.0001` 的浮点波段用于展示和 Sen 斜率。

```javascript
function prepNDVI(image) {
  var ndvi = image.select('NDVI');
  var qa = image.select('SummaryQA');
  var mask = qa.lte(1)
    .and(ndvi.gte(-2000))
    .and(ndvi.lte(10000));
  return ndvi.updateMask(mask)
    .rename('NDVI_raw')
    .clipToCollection(regionFC)
    .setDefaultProjection(nativeProj)
    .copyProperties(image, ['system:time_start']);
}
```

## 3. GPP 改造口径（山西案例）

- 产品：`MODIS/061/MOD17A2HGF`，500 m，8 天，gap-filled 累计 GPP；
- 波段：`Gpp`（另有 `PsnNet`、`Psn_QC`）；
- 时段值：原始整数 ∈ `[0, 3000]` 视为有效，保留产品掩膜，**缺测不补 0**；
- 年度值：日历年 1 月 1 日—12 月 31 日，**46 个 8 天时段全部有效**才累加；
- 物理值：原始整数 × 0.0001 = kg C/m²/时段；全年累加 × 0.1 = **g C/m²/year**；
- 时间：2001—2022，共 22 年；
- 存储：年度累加值保持 Int32，分析时再乘 0.1。

```javascript
var source = ee.ImageCollection('MODIS/061/MOD17A2HGF')
  .filterDate(start, end)
  .filterBounds(roi);

function prepGPP(image) {
  var raw = image.select('Gpp');
  return raw.updateMask(raw.gte(0).and(raw.lte(3000)))
    .rename('GPP_raw')
    .setDefaultProjection(nativeProj)
    .copyProperties(image, ['system:time_start']);
}

var annual = ee.ImageCollection.fromImages(
  years.map(function (y) {
    y = ee.Number(y);
    var start = ee.Date.fromYMD(y, 1, 1);
    var end = start.advance(1, 'year');
    var imgs = source.filterDate(start, end).map(prepGPP);
    var rawAnnual = imgs.sum().rename('GPP_raw_annual');
    var complete = imgs.count().eq(46);
    return rawAnnual.updateMask(complete)
      .set('year', y)
      .set('system:time_start', start.millis());
  })
);
```

`imgs.count().eq(46)` 是关键：某像元只要缺一个时段，整年作废，不参与趋势计算。
这与“缺测补 0”完全不同——补 0 会制造虚假的下降趋势。

### 为什么不用 MOD17A3HGF 年产品？

`MODIS/061/MOD17A3HGF` 确实是年 GPP 产品，也包含 `Gpp` 波段，可以作为替代方案。
本项目仍选 8 天产品，原因是：

- 年度累加的定义、QA 和缺测处理完全由我们控制，口径可写清楚；
- 与 NDVI 原版一样保持“从原始周期产品合成年度值”的流程，便于方法对照；
- 8 天产品可以同时做生长季/非生长季拆分，后续扩展更灵活。

如果直接用年产品，必须在方法节写明产品版本、QC 策略和比例因子，不能混用两套口径。

## 4. 完整序列掩膜

无论 NDVI 还是 GPP，本仓库都要求**每个像元在全部年份有效**：

```javascript
var count = annual.select('value').count();
var completeMask = count.eq(N_YEARS).selfMask();
```

原因：

- Sen 斜率要求年份对可比；
- MK 的 $n$ 必须一致，否则 $S$ 与方差不可比；
- 插值/补 0 会引入人为趋势。

代价是覆盖率下降。山西案例中完整序列像元占矩形网格的 57.9%，
但对应省内有效面积 151,728.7 km²，占全省 96.8%——损失主要在省界外的矩形角。

## 5. 不同年度指标会得到不同结论

同一份 8 天 GPP，若改用“生长季最大值”，趋势会明显更平（峰值受单季限制），
显著像元比例也会下降。这不是算法问题，而是科学问题变了。

因此论文必须写清楚四件事：

1. 产品名称、版本、空间/时间分辨率；
2. 年度指标的合成方式（max / mean / sum / median）；
3. 时间窗口（日历年 or 生长季；起止月日）；
4. 缺测与 QA 策略（完整序列 / 允许缺测 / 插值）。

缺任何一项，结果都不可复现。

## 6. 其他可选指标

| 指标 | 适用 | 注意 |
|---|---|---|
| 生长季最大值 | 植被长势峰值 | 对单期噪声敏感，需 QA |
| 生长季均值 | 整体长势 | 受生长季长度定义影响 |
| 年累加 | GPP、NPP、ET、降水 | 缺测不能补 0 |
| 年 median | 抗异常值 | 适合有极端值的状态量 |
| 分位数（p90） | 极端事件/峰值 | 需足够多的年内观测 |
| 物候指标（SOS/EOS） | 生长季起止 | 与趋势方法分开讨论 |

本仓库的代码把年度指标留在 GEE 侧（`01_check_annual_gpp.js` / `02_export...js`），
本地脚本只接收已经合成好的年度栅格，避免“年度口径”和“趋势算法”耦合。
