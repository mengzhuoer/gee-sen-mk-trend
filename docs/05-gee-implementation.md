# 05 · GEE 实现：从年度集合到趋势图

本仓库的 GEE 代码分两类：

- `code/gee/sen_mk_ndvi.js`：NDVI 完整流程（对应原教程方法，2001—2024，MOD13A1）；
- `code/gee/gpp_shanxi/`：山西 GPP 案例的**分阶段**脚本（01—05，已在 GEE 实跑通过）。

下面按计算顺序说明关键实现与坑位。

## 1. 研究区

```javascript
var regionFC = ee.FeatureCollection('projects/xxx/assets/your_roi');
var roi = regionFC.geometry();
```

若用 GAUL 行政区，**不要按拼音精确匹配省名**。GAUL 的山西写作 `Shanxi Sheng`，
按 `ADM1_NAME == 'Shanxi'` 会得到空集合，后续所有 `first()` 都变成 null。
稳妥做法是用省会坐标空间选取，并打印实际名称：

```javascript
var taiyuan = ee.Geometry.Point([112.549, 37.871]);
var regionFC = ee.FeatureCollection('FAO/GAUL/2015/level1')
  .filterBounds(taiyuan);
print('要素数', regionFC.size());
print('名称', regionFC.aggregate_array('ADM1_NAME'));
```

## 2. 投影：不要让常量影像落到 1° 默认网格

GEE 的 `ee.Image.constant()` 默认是 1° 的 WGS84 网格。一旦它和 500 m 影像参与
同一个 `addBands` / `where` / 数组运算，整个计算可能被重采样到 1°，结果直接错。

两种防护：

1. 所有派生影像显式 `setDefaultProjection(nativeProj)`；
2. 常量波段从已有影像派生，继承其网格：

```javascript
var yearBand = image.select('value').multiply(0).add(year);
```

文章代码大量使用 `setDefaultProjection(nativeProj)`，不是冗余，是必需的。

## 3. 掩膜：`unmask(0)` 是趋势分析的头号杀手

```javascript
// 错误：把研究区外的缺测变成 0，参与统计后制造虚假趋势
image.unmask(0)

// 正确：保持掩膜，最后用完整序列掩膜筛像元
image.updateMask(completeMask)
```

原文在 MK 的 sign 影像处特意注明“no unmask(0)”，原因就在此。
GPP 的 0 是合法观测值，缺测是掩膜，两者必须区分。

## 4. Sen 斜率

```javascript
var senInput = annual.map(function (img) {
  var y = ee.Number(img.get('year'));
  var yearBand = img.select('value').multiply(0).add(y).rename('year');
  return yearBand.addBands(img.select('value'))
    .setDefaultProjection(nativeProj);
});

var sen = senInput.select(['year', 'value'])
  .reduce(ee.Reducer.sensSlope());

var senSlope = sen.select('slope').rename('SenSlope')
  .updateMask(completeMask);
```

`ee.Reducer.sensSlope()` 输出 `slope` 和 `offset`；我们只用 `slope`。
斜率单位是“value 单位/年”。

## 5. MK 的 S：Join.saveAll

```javascript
var afterFilter = ee.Filter.lessThan({
  leftField: 'system:time_start',
  rightField: 'system:time_start'
});

var joined = ee.ImageCollection(
  ee.Join.saveAll('after').apply({
    primary: mkCollection,
    secondary: mkCollection,
    condition: afterFilter
  })
);

function sign(i, j) {
  i = ee.Image(i);
  j = ee.Image(j);
  return j.neq(i)
    .multiply(j.subtract(i).clamp(-1, 1))
    .rename('sign')
    .toInt();
}

var signImages = ee.ImageCollection(
  joined.map(function (current) {
    var after = ee.ImageCollection.fromImages(current.get('after'));
    return after.map(function (next) {
      return sign(current, next);
    });
  }).flatten()
);

var S = signImages.reduce(ee.Reducer.sum(), 2).rename('S');
```

`j.neq(i)` 处理 ties（相等记 0）；`clamp(-1,1)` 把差值变成符号。
`reduce(ee.Reducer.sum(), 2)` 的第二个参数是并行尺度，不是“按第 2 维求和”。

## 6. ties 修正的数组写法

```javascript
var groups = mkCollection.map(function (i) {
  i = ee.Image(i);
  var matches = mkCollection.map(function (j) {
    return i.eq(ee.Image(j));
  }).sum();
  // 非重复年份 -> 0（哨兵组），重复年份保留原值
  return i.multiply(matches.gt(1));
});

var n = mkCollection.count();
var varS = tieFactor(n)
  .subtract(groupFactorSum)
  .divide(18)
  .rename('VarS');
```

`groupSizesFromArray` 排序后统计连续相同值的游程长度，并丢弃第一组（哨兵 0 组）。
真实 0 值会被并入哨兵组，这是该写法的已知偏差；本地精确版见
[02-mann-kendall.md](02-mann-kendall.md) 第 6 节。

## 7. Z 与 p

```javascript
var Z = S.multiply(0)
  .where(S.gt(0), S.subtract(1).divide(varS.sqrt()))
  .where(S.lt(0), S.add(1).divide(varS.sqrt()))
  .updateMask(varS.gt(0))
  .updateMask(completeMask);

function normalCDF(z) {
  return ee.Image(z).divide(Math.sqrt(2)).erf().add(1).multiply(0.5);
}

var pValue = ee.Image(1)
  .subtract(normalCDF(Z.abs()))
  .multiply(2)
  .clamp(0, 1);
```

`VarS=0` 的常量像元会被 `varS.gt(0)` 掩掉；若希望显式归为稳定类，
可在本地实现中定义 $Z=0,p=1$（本仓库采用后者）。

## 8. 容量限制与分阶段

山西案例的单体脚本在同时请求

- 22 年逐年累加；
- 231 个年份对的 MK；
- 全省分位数、年际曲线、面积统计；
- 图层渲染

时触发：

```
Earth Engine capacity exceeded
```

处理策略：

1. **拆阶段**：检查（01）→ 导出/下载（02、03、05）→ 指纹（04）→ 本地趋势（06）；
2. **每年一个独立小计算**，避免共同上游被重复求值；
3. **趋势计算搬到本地**，7—10 秒完成，且完全可复现；
4. 若必须留在 GEE，可先 `Export` 年度集合，再对资产做 Sen/MK，减少一次请求的依赖链。

## 9. 导出与下载

### Export 到 Drive / Assets

```javascript
Export.image.toDrive({
  image: out,
  description: 'Shanxi_GPP_SenMK_2001_2022',
  folder: 'GEE_SenMK',
  fileNamePrefix: 'Shanxi_GPP_SenMK_2001_2022',
  region: roi,
  scale: 500,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
```

大区域建议 `crs: 'EPSG:4326'` + 固定 `crsTransform`，保证所有年份同一网格；
否则不同年份可能因默认投影/分块产生细微网格差异，叠栈时报错或错位。

### getDownloadURL

```javascript
var url = image.getDownloadURL({
  region: roi,
  scale: 500,
  format: 'GEO_TIFF'
});
print(url);
```

已知坑：

- `filename` 参数会被忽略，下载文件名是哈希；
- `crs_transform` 不接受 ee 对象，必须是普通 JS 数组/字符串；
- 链接是签名 URL，生成后尽快下载，过期重跑；
- 单次请求有大小/像元上限，大区域会失败。

文件变成哈希名后，用“多点指纹”认领年份，见
[07-pitfalls.md](07-pitfalls.md) 与独立仓库 `yaogan-zhiwen`。

## 10. 运行前检查清单

- [ ] `print` 研究区要素数 = 1，名称正确；
- [ ] `print` 影像数 = 年数 × 每年时段数；
- [ ] `print` 波段名、比例因子与 catalog 一致；
- [ ] 每年时段数齐全（NDVI 23 期、GPP 46 期，按产品实际值核对）；
- [ ] `print` 完整序列像元数 / 覆盖率；
- [ ] 先跑单年，确认导出/下载网格，再批量；
- [ ] 同一批文件 CRS + transform 完全一致；
- [ ] 结果图与年际曲线同时检查，避免“图对但数错”。
