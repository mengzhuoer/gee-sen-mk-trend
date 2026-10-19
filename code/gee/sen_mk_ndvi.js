/***************************************************************
 * Sen's slope + Mann–Kendall 长时序栅格趋势分析（NDVI 示例）
 *
 * 方法对应文档：
 *   docs/01-sen-slope.md
 *   docs/02-mann-kendall.md
 *   docs/03-trend-classification.md
 *   docs/04-annual-metric.md
 *   docs/05-gee-implementation.md
 *
 * 数据：MODIS/061/MOD13A1，500 m，16 天
 * 年度指标：每年 4—10 月最大 NDVI（演示口径，可按科学问题替换）
 * 分级：1 显著下降 / 2 不显著下降 / 3 基本稳定 / 4 不显著上升 / 5 显著上升
 *
 * 说明：
 * - 本脚本是方法复刻与整理，不转载原教程正文；
 * - 大区域若报 `Earth Engine capacity exceeded`，请拆阶段或改用本地计算
 *   （见 code/python/sen_mk_local.py）；
 * - 运行前把 ROI_ASSET 换成你自己的矢量资产，或把 ROI_ASSET 设为 null
 *   使用 GAUL 行政区（按省会坐标空间选取）。
 ***************************************************************/

// ==================== 1. 参数 ====================

var START_YEAR = 2001;
var END_YEAR = 2024;
var N_YEARS = END_YEAR - START_YEAR + 1;
var SCALE = 500;
var ALPHA = 0.05;
var EXPORT_RESULTS = false;

// 替换为你的研究区资产；设为 null 时用 GAUL 省级行政区（太原点 → Shanxi Sheng）。
var ROI_ASSET = null;
var GAUL_POINT = [112.549, 37.871];

var regionFC = ROI_ASSET
  ? ee.FeatureCollection(ROI_ASSET)
  : ee.FeatureCollection('FAO/GAUL/2015/level1')
      .filterBounds(ee.Geometry.Point(GAUL_POINT));

var roi = regionFC.geometry();
var years = ee.List.sequence(START_YEAR, END_YEAR);

print('研究区要素数（应为 1）：', regionFC.size());
print('研究区名称：', regionFC.aggregate_array('ADM1_NAME'));

// ==================== 2. 数据与 QA ====================

var mod13 = ee.ImageCollection('MODIS/061/MOD13A1')
  .filterDate(START_YEAR + '-01-01', (END_YEAR + 1) + '-01-01')
  .filterBounds(roi);

print('MOD13A1 影像总数：', mod13.size());

// 原生投影：所有派生影像都显式设置，避免落到 1° 默认网格。
var nativeProj = ee.Image(mod13.first()).select('NDVI').projection();
print('原生投影：', nativeProj);

// SummaryQA: 0 = good, 1 = marginal, 2 = snow/ice, 3 = cloudy.
function prepNDVI(image) {
  image = ee.Image(image);
  var ndvi = image.select('NDVI');
  var qa = image.select('SummaryQA');

  var mask = qa.lte(1)
    .and(ndvi.gte(-2000))
    .and(ndvi.lte(10000));

  return ndvi
    .updateMask(mask)
    .rename('NDVI_raw')
    .clipToCollection(regionFC)
    .setDefaultProjection(nativeProj)
    .copyProperties(image, ['system:time_start']);
}

// ==================== 3. 年度 NDVI ====================

var annualNDVI = ee.ImageCollection.fromImages(
  years.map(function (y) {
    y = ee.Number(y);

    var start = ee.Date.fromYMD(y, 4, 1);
    var end = ee.Date.fromYMD(y, 11, 1);

    var raw = mod13
      .filterDate(start, end)
      .map(prepNDVI)
      .max()
      .rename('NDVI_raw')
      .toInt()
      .setDefaultProjection(nativeProj)
      .clipToCollection(regionFC);

    var value = raw
      .multiply(0.0001)
      .rename('NDVI')
      .toFloat()
      .setDefaultProjection(nativeProj)
      .clipToCollection(regionFC);

    return value.addBands(raw)
      .setDefaultProjection(nativeProj)
      .clipToCollection(regionFC)
      .set({
        year: y,
        'system:time_start': ee.Date.fromYMD(y, 7, 1).millis()
      });
  })
);

print('年度 NDVI 影像数：', annualNDVI.size());

// 完整序列掩膜：N_YEARS 年全部有效才参与分析。
var count = annualNDVI.select('NDVI').count()
  .setDefaultProjection(nativeProj)
  .clipToCollection(regionFC);

var completeMask = count.eq(N_YEARS)
  .selfMask()
  .setDefaultProjection(nativeProj)
  .clipToCollection(regionFC);

// ==================== 4. Sen's slope ====================

var senInput = annualNDVI.map(function (image) {
  image = ee.Image(image);
  var y = ee.Number(image.get('year'));

  // 从 NDVI 波段派生年份波段，继承同一网格与掩膜。
  var yearBand = image.select('NDVI')
    .multiply(0)
    .add(y)
    .rename('year')
    .toFloat();

  return yearBand
    .addBands(image.select('NDVI'))
    .setDefaultProjection(nativeProj)
    .clipToCollection(regionFC)
    .copyProperties(image, ['system:time_start']);
});

var sen = senInput
  .select(['year', 'NDVI'])
  .reduce(ee.Reducer.sensSlope());

var senSlope = sen.select('slope')
  .rename('SenSlope')
  .updateMask(completeMask)
  .clipToCollection(regionFC);

// ==================== 5. MK 的 S ====================

var mkCollection = annualNDVI.map(function (image) {
  return ee.Image(image)
    .select('NDVI_raw')
    .toInt()
    .setDefaultProjection(nativeProj)
    .clipToCollection(regionFC)
    .copyProperties(image, ['system:time_start']);
});

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
    .toInt()
    .setDefaultProjection(nativeProj)
    .clipToCollection(regionFC);
}

// 注意：不要 unmask(0)。完整序列像元之外保持掩膜，
// 否则会被填 0 并重新出现在结果里。
var signImages = ee.ImageCollection(
  joined.map(function (current) {
    current = ee.Image(current);
    var after = ee.ImageCollection.fromImages(current.get('after'));
    return after.map(function (next) {
      return sign(current, next);
    });
  }).flatten()
);

var S = signImages
  .reduce(ee.Reducer.sum(), 2)
  .rename('S')
  .setDefaultProjection(nativeProj)
  .updateMask(completeMask)
  .clipToCollection(regionFC);

// ==================== 6. ties 修正方差 ====================

// 该年取值出现次数 > 1 的保留原值，否则置 0（哨兵组）。
var groups = mkCollection.map(function (i) {
  i = ee.Image(i);
  var matches = mkCollection.map(function (j) {
    return i.eq(ee.Image(j));
  }).sum();

  return i.multiply(matches.gt(1))
    .setDefaultProjection(nativeProj)
    .clipToCollection(regionFC);
});

function groupSizesFromArray(array) {
  var length = array.arrayLength(0);

  // 原生网格上的常量 1，避免 1° 默认网格。
  var one = ee.Image.constant(1)
    .setDefaultProjection(nativeProj)
    .clipToCollection(regionFC);

  var indices = one
    .arrayRepeat(0, length)
    .arrayAccum(0, ee.Reducer.sum())
    .toArray(1);

  var sorted = array.arraySort();
  var left = sorted.arraySlice(0, 1);
  var right = sorted.arraySlice(0, 0, -1);

  // 最后一个游程没有“右邻居”，显式补一个结束标记。
  var last = ee.Image.constant(ee.Array([[1]]))
    .setDefaultProjection(nativeProj)
    .clipToCollection(regionFC);

  var endMask = left.neq(right).arrayCat(last, 0);
  var runIndices = indices.arrayMask(endMask);

  // 相邻游程结束位置之差 = 游程长度；
  // 第一个游程是 0 哨兵组，恰好被丢掉。
  return runIndices.arraySlice(0, 1)
    .subtract(runIndices.arraySlice(0, 0, -1));
}

function tieFactor(image) {
  return image.expression('b() * (b() - 1) * (2 * b() + 5)');
}

var groupSizes = groupSizesFromArray(groups.toArray());
var groupFactorSum = tieFactor(groupSizes)
  .arrayReduce(ee.Reducer.sum(), [0])
  .arrayGet([0, 0]);

var n = mkCollection.count()
  .setDefaultProjection(nativeProj)
  .clipToCollection(regionFC);

var varS = tieFactor(n)
  .subtract(groupFactorSum)
  .divide(18)
  .rename('VarS')
  .toFloat()
  .updateMask(completeMask)
  .clipToCollection(regionFC);

// ==================== 7. Z 与双侧 p ====================

var Z = S.multiply(0)
  .where(S.gt(0), S.subtract(1).divide(varS.sqrt()))
  .where(S.lt(0), S.add(1).divide(varS.sqrt()))
  .rename('Z')
  .updateMask(varS.gt(0))
  .updateMask(completeMask)
  .clipToCollection(regionFC);

function normalCDF(z) {
  return ee.Image(z)
    .divide(Math.sqrt(2))
    .erf()
    .add(1)
    .multiply(0.5);
}

var pValue = ee.Image(1)
  .subtract(normalCDF(Z.abs()))
  .multiply(2)
  .clamp(0, 1)
  .rename('P')
  .updateMask(completeMask)
  .clipToCollection(regionFC);

// ==================== 8. 五级趋势 ====================

var trendClass = senSlope.multiply(0).toByte()
  .where(senSlope.lt(0).and(pValue.lt(ALPHA)), 1)
  .where(senSlope.lt(0).and(pValue.gte(ALPHA)), 2)
  .where(senSlope.eq(0), 3)
  .where(senSlope.gt(0).and(pValue.gte(ALPHA)), 4)
  .where(senSlope.gt(0).and(pValue.lt(ALPHA)), 5)
  .rename('TrendClass')
  .updateMask(completeMask)
  .clipToCollection(regionFC);

// ==================== 9. 地图与图例 ====================

var trendPalette = ['A50026', 'F46D43', 'D9D9D9', '66BD63', '006837'];

Map.centerObject(regionFC, 6);

Map.addLayer(senSlope,
  {min: -0.01, max: 0.01,
   palette: ['8E0152', 'DE77AE', 'F7F7F7', '7FBC41', '276419']},
  'NDVI - Sen slope', false);

Map.addLayer(Z,
  {min: -3, max: 3,
   palette: ['B2182B', 'EF8A62', 'F7F7F7', '67A9CF', '2166AC']},
  'NDVI - MK Z', false);

Map.addLayer(pValue,
  {min: 0, max: 0.1, palette: ['006837', 'FFFFBF', 'A50026']},
  'NDVI - two-sided p', false);

Map.addLayer(trendClass,
  {min: 1, max: 5, palette: trendPalette},
  'NDVI - Sen+MK trend class', true);

Map.addLayer(
  regionFC.style({color: '222222', fillColor: '00000000', width: 2}),
  {},
  'ROI boundary', true);

var legend = ui.Panel({style: {position: 'bottom-left', padding: '8px 12px'}});
legend.add(ui.Label({
  value: 'Sen + MK 趋势分级',
  style: {fontWeight: 'bold', fontSize: '14px', margin: '0 0 6px 0'}
}));

var legendText = [
  '显著退化（Slope < 0，p < 0.05）',
  '不显著退化（Slope < 0，p ≥ 0.05）',
  '基本稳定（Slope = 0）',
  '不显著改善（Slope > 0，p ≥ 0.05）',
  '显著改善（Slope > 0，p < 0.05）'
];

for (var k = 0; k < 5; k++) {
  legend.add(ui.Panel({
    widgets: [
      ui.Label({style: {
        backgroundColor: '#' + trendPalette[k],
        padding: '8px',
        margin: '0 6px 4px 0'
      }}),
      ui.Label({value: legendText[k], style: {margin: '0 0 4px 0'}})
    ],
    layout: ui.Panel.Layout.Flow('horizontal')
  }));
}
Map.add(legend);

// ==================== 10. 年际曲线 ====================

var seriesFC = ee.FeatureCollection(
  years.map(function (y) {
    y = ee.Number(y);

    var image = ee.Image(
      annualNDVI.filter(ee.Filter.eq('year', y)).first()
    );

    var mean = image.select('NDVI').reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: roi,
      scale: SCALE,
      maxPixels: 1e13,
      tileScale: 4
    }).get('NDVI');

    return ee.Feature(null, {year: y, NDVI: mean});
  })
);

var seriesChart = ui.Chart.feature.byFeature(seriesFC, 'year', ['NDVI'])
  .setChartType('LineChart')
  .setOptions({
    title: START_YEAR + '–' + END_YEAR + ' 年研究区最大 NDVI 年际变化',
    hAxis: {title: 'Year'},
    vAxis: {title: 'NDVI'},
    lineWidth: 2,
    pointSize: 4,
    legend: {position: 'none'}
  });
print(seriesChart);

// ==================== 11. 面积统计 ====================

var classNames = ee.Dictionary({
  '1': '显著退化',
  '2': '不显著退化',
  '3': '基本稳定',
  '4': '不显著改善',
  '5': '显著改善'
});

var areaImage = ee.Image.pixelArea().divide(1e6).rename('area_km2')
  .addBands(trendClass.rename('class'))
  .clipToCollection(regionFC);

var grouped = areaImage.reduceRegion({
  reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'class'}),
  geometry: roi,
  scale: SCALE,
  maxPixels: 1e13,
  tileScale: 4
});

var groupsList = ee.List(grouped.get('groups'));
var areaFC = ee.FeatureCollection(groupsList.map(function (item) {
  item = ee.Dictionary(item);
  var cls = ee.Number(item.get('class')).int();
  return ee.Feature(null, {
    class: cls,
    trend: classNames.get(cls.format()),
    area_km2: item.get('sum')
  });
})).sort('class');

var totalArea = ee.Number(areaFC.aggregate_sum('area_km2'));
areaFC = areaFC.map(function (f) {
  return f.set('percentage',
    ee.Number(f.get('area_km2')).divide(totalArea).multiply(100));
});

print('NDVI 趋势类型面积：', areaFC);

var areaChart = ui.Chart.feature.byFeature(areaFC, 'trend', ['percentage'])
  .setChartType('ColumnChart')
  .setOptions({
    title: 'NDVI 趋势类型面积占比',
    hAxis: {title: '趋势类型'},
    vAxis: {title: 'Percentage (%)'},
    legend: {position: 'none'}
  });
print(areaChart);

// ==================== 12. 导出 ====================

if (EXPORT_RESULTS) {
  var out = senSlope
    .addBands(Z)
    .addBands(pValue)
    .addBands(trendClass)
    .clipToCollection(regionFC);

  Export.image.toDrive({
    image: out,
    description: 'NDVI_SenMK_' + START_YEAR + '_' + END_YEAR,
    folder: 'GEE_SenMK',
    fileNamePrefix: 'NDVI_SenMK_' + START_YEAR + '_' + END_YEAR,
    region: roi,
    scale: SCALE,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  });
}

print('完成：' + START_YEAR + '–' + END_YEAR + ' NDVI Sen + MK');
