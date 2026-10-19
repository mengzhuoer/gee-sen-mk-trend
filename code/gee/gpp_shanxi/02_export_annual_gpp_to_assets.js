/**************************************************************
 * 第二步：将山西省年度 GPP 导出到 GEE Assets
 *
 * 数据：MODIS/061/MOD17A2HGF
 * 指标：全年各时段 GPP 累加
 * 当前边界：GAUL 山西省
 *
 * 首次只导出 2022 年，成功后再分批导出其余年份。
 *
 * 导出波段：GPP_raw_annual，Int32
 * 物理单位转换：
 * GPP_raw_annual × 0.1 = g C/m²/year
 *
 * 不进行 Sen、MK、图表或省域统计。
 **************************************************************/


// ==================== 1. 参数 ====================

// 首次只导出 2022 年。
// 验证成功后，可改成小批次，例如 2001—2005。
var EXPORT_START_YEAR = 2022;
var EXPORT_END_YEAR = 2022;

// 当前云项目。
// 如导出权限报错，应检查此项目的资产写入权限。
var ASSET_ROOT = 'projects/songyuan/assets';

// 明确标注当前使用 GAUL 边界。
var ASSET_PREFIX = 'Shanxi_GAUL_GPP_raw_annual_';


// ==================== 2. 山西省边界 ====================

// 用太原坐标选择省级边界。
// 已验证对应 GAUL 的 Shanxi Sheng。
var taiyuanPoint = ee.Geometry.Point([112.549, 37.871]);

var regionFC = ee.FeatureCollection(
  'FAO/GAUL/2015/level1'
).filterBounds(taiyuanPoint);

var roi = regionFC.geometry();

Map.clear();
Map.centerObject(taiyuanPoint, 7);

print('第二步：年度 GPP 导出到 Assets。');
print('本次导出年份：',
  EXPORT_START_YEAR + '—' + EXPORT_END_YEAR);

print('研究区要素数，预期为 1：',
  regionFC.size());

print('省级名称，预期为 Shanxi Sheng：',
  regionFC.aggregate_array('ADM1_NAME'));


// ==================== 3. 读取数据 ====================

var source = ee.ImageCollection(
  'MODIS/061/MOD17A2HGF'
)
  .filterDate(
    ee.Date.fromYMD(EXPORT_START_YEAR, 1, 1),
    ee.Date.fromYMD(EXPORT_END_YEAR + 1, 1, 1)
  )
  .filterBounds(roi)
  .sort('system:time_start');

print('本批次 8 天影像数量：', source.size());

var nativeProj = ee.Image(source.first())
  .select('Gpp')
  .projection();

// 只获取投影元数据，不计算全省栅格统计。
// 后续导出使用相同的原生网格。
var projectionInfo = nativeProj.getInfo();

print('原生投影：', nativeProj);


// ==================== 4. 有效值处理 ====================

// 保留产品已有掩膜。
// 真实 0 值保留；范围外填充值排除。
// 不将缺失时段填成 0。
function prepareGPP(image) {
  image = ee.Image(image);

  var raw = image.select('Gpp');

  return raw
    .updateMask(raw.gte(0).and(raw.lte(3000)))
    .rename('raw')
    .toInt32()
    .copyProperties(image, ['system:time_start']);
}


// ==================== 5. 构建单年累计 GPP ====================

function buildAnnualRaw(year) {
  var start = ee.Date.fromYMD(year, 1, 1);
  var end = start.advance(1, 'year');

  var yearSource = source.filterDate(start, end);
  var prepared = yearSource.map(prepareGPP);

  var nImages = yearSource.size();

  var validCount = prepared.select('raw')
    .count()
    .setDefaultProjection(nativeProj);

  // 本年必须有完整的 46 个时段。
  var yearIsComplete = ee.Image.constant(
    nImages.eq(46)
  ).setDefaultProjection(nativeProj);

  // 每个像元必须在全部时段都有有效值。
  var annualMask = validCount.eq(nImages)
    .and(yearIsComplete);

  var annualRaw = prepared.select('raw')
    .sum()
    .rename('GPP_raw_annual')
    .toInt32()
    .setDefaultProjection(nativeProj)
    .updateMask(annualMask)
    .clipToCollection(regionFC);

  return annualRaw.set({
    year: year,
    n_images: nImages,
    source_product: 'MODIS/061/MOD17A2HGF',
    annual_metric: 'calendar-year sum of interval GPP',
    boundary_source: 'FAO/GAUL/2015/level1',
    boundary_selection: 'province containing Taiyuan',
    value_description: 'annual sum of original Gpp integers',
    to_gC_m2_year_factor: 0.1,
    completeness_rule: '46 intervals, all valid per pixel',
    'system:time_start': start.millis()
  });
}


// ==================== 6. 创建导出任务 ====================

// JavaScript 客户端循环：每年创建一个独立导出任务。
// 点击编辑器 Run 不会自动提交这些任务。
for (
  var year = EXPORT_START_YEAR;
  year <= EXPORT_END_YEAR;
  year++
) {
  var annualRaw = buildAnnualRaw(year);

  var assetName = ASSET_PREFIX + year;
  var assetId = ASSET_ROOT + '/' + assetName;

  print(year + ' 年时段数量，预期为 46：',
    annualRaw.get('n_images'));

  print('目标资产：', assetId);

  Export.image.toAsset({
    image: annualRaw,
    description: assetName,
    assetId: assetId,
    region: roi,

    // 保持 MODIS 原生网格，不另外指定 scale。
    crs: projectionInfo.crs,
    crsTransform: projectionInfo.transform,

    // 金字塔显示采用平均值。
    // 后续趋势计算仍应使用原生分辨率。
    pyramidingPolicy: {
      '.default': 'mean'
    },

    maxPixels: 1e13
  });
}


// ==================== 7. 显示边界 ====================

// 本步骤不添加年度 GPP 地图，
// 避免预览与导出重复请求年度栅格计算。
Map.addLayer(
  regionFC.style({
    color: '222222',
    fillColor: '00000000',
    width: 2
  }),
  {},
  '山西省 GAUL 边界',
  true
);

print(
  '任务已创建。请到 Tasks 面板找到本次任务，' +
  '点击 RUN 并确认提交。'
);

print(
  '只有任务状态变为 COMPLETED，才表示导出成功。' +
  '此 Console 提示不代表导出已经完成。'
);