/**************************************************************
 * 第一步：山西省边界、数据完整性与年度 GPP 检查
 *
 * 时间：2001—2022 年
 * 数据：MODIS/061/MOD17A2HGF
 * 指标：8 天累计 GPP 按全年求和
 * 年度单位：g C/m²/year
 *
 * 本步骤不运行：
 * Sen、MK、全省统计、图表、导出
 **************************************************************/


// ==================== 1. 参数 ====================

var START_YEAR = 2001;
var END_YEAR = 2022;
var N_YEARS = END_YEAR - START_YEAR + 1;

// 本次只显示这一年的地图。
var DISPLAY_YEAR = 2022;

// 当前先用 GAUL 验证流程。
// 后面使用 2023 版边界时，再填写真实资产 ID。
var USE_UPLOADED_ASSET = false;

var ASSET_PATH =
  'projects/songyuan/assets/填入真实的山西省资产名';


// ==================== 2. 山西省边界 ====================

// 太原市中心，用空间位置选择省级行政区。
// 已在你的账号中验证：选中的是 Shanxi Sheng。
var taiyuanPoint = ee.Geometry.Point([112.549, 37.871]);

var gaulShanxi = ee.FeatureCollection(
  'FAO/GAUL/2015/level1'
).filterBounds(taiyuanPoint);

var regionFC = USE_UPLOADED_ASSET
  ? ee.FeatureCollection(ASSET_PATH)
  : gaulShanxi;

var roi = regionFC.geometry();

Map.clear();
Map.centerObject(taiyuanPoint, 7);

print('第一步：检查年度 GPP，尚未进行趋势分析。');

print('研究区要素数，预期为 1：',
  regionFC.size());

print('实际省级名称，GAUL 预期为 Shanxi Sheng：',
  regionFC.aggregate_array('ADM1_NAME'));


// ==================== 3. 原始数据 ====================

var source = ee.ImageCollection(
  'MODIS/061/MOD17A2HGF'
)
  .filterDate(
    ee.Date.fromYMD(START_YEAR, 1, 1),
    ee.Date.fromYMD(END_YEAR + 1, 1, 1)
  )
  .filterBounds(roi)
  .sort('system:time_start');

print('8 天影像总数，预期为 1012：',
  source.size());

print('首张影像波段，预期包含 Gpp：',
  ee.Image(source.first()).bandNames());

var nativeProj = ee.Image(source.first())
  .select('Gpp')
  .projection();


// ==================== 4. 年份完整性检查 ====================

var years = ee.List.sequence(START_YEAR, END_YEAR);

var yearCheck = ee.FeatureCollection(
  years.map(function(year) {
    year = ee.Number(year);

    var start = ee.Date.fromYMD(year, 1, 1);
    var end = start.advance(1, 'year');

    return ee.Feature(null, {
      year: year,
      n_images: source.filterDate(start, end).size()
    });
  })
).sort('year');

print('逐年影像数量，请展开检查：',
  yearCheck);

print('时段数量不是 46 的年份，预期为空：',
  yearCheck.filter(
    ee.Filter.neq('n_images', 46)
  ));


// ==================== 5. 原始 GPP 有效值处理 ====================

// 产品已有掩膜继续保留。
// 有效原始值范围为 0—3000。
// 真实 0 值是合法观测，不当作缺失值。
//
// 本步骤不额外施加 Psn_QC 位筛选。
// 正式研究时需要明确质量控制口径。
function prepareGPP(image) {
  image = ee.Image(image);

  var raw = image.select('Gpp');

  return raw
    .updateMask(raw.gte(0).and(raw.lte(3000)))
    .rename('raw')
    .toInt32()
    .copyProperties(image, ['system:time_start']);
}


// ==================== 6. 构建年度 GPP ====================

// GEE 按需计算。
// 构建 22 年集合不等于立即计算全部年度栅格；
// 本次地图只请求 DISPLAY_YEAR 的结果。
var annualGPP = ee.ImageCollection.fromImages(
  years.map(function(year) {
    year = ee.Number(year);

    var start = ee.Date.fromYMD(year, 1, 1);
    var end = start.advance(1, 'year');

    var yearSource = source.filterDate(start, end);
    var prepared = yearSource.map(prepareGPP);

    var nImages = yearSource.size();

    // 像元在该年具有有效值的时段数量。
    var validCount = prepared.select('raw')
      .count()
      .setDefaultProjection(nativeProj);

    // 一年正常有 46 个时段。
    // 最后一个时段不足 8 天，仍直接使用产品累计量，
    // 不额外乘以 8。
    var yearIsComplete = ee.Image.constant(
      nImages.eq(46)
    ).setDefaultProjection(nativeProj);

    // 严格要求：全年时段完整，且像元每个时段都有值。
    // 不把缺测时段填为 0。
    var annualMask = validCount.eq(nImages)
      .and(yearIsComplete);

    // 全年累计原始整数，后续供 MK 判断重复值。
    var rawAnnual = prepared.select('raw')
      .sum()
      .rename('GPP_raw_annual')
      .toInt32()
      .setDefaultProjection(nativeProj)
      .updateMask(annualMask);

    // raw × 0.0001 kg C/m² × 1000
    // = raw × 0.1 g C/m²
    var annualValue = rawAnnual
      .multiply(0.1)
      .rename('GPP')
      .toFloat();

    return annualValue
      .addBands(rawAnnual)
      .clipToCollection(regionFC)
      .set({
        year: year,
        n_images: nImages,
        'system:time_start': start.millis()
      });
  })
).sort('system:time_start');

print('年度集合数量，预期为 22：',
  annualGPP.size());


// ==================== 7. 只显示一个年度 ====================

var displayImage = ee.Image(
  annualGPP
    .filter(ee.Filter.eq('year', DISPLAY_YEAR))
    .first()
);

Map.addLayer(
  displayImage.select('GPP'),
  {
    min: 0,
    max: 2000,
    palette: [
      'FFFFCC',
      'C2E699',
      '78C679',
      '31A354',
      '006837'
    ]
  },
  DISPLAY_YEAR + ' 年度 GPP，g C/m²/year',
  true
);

Map.addLayer(
  regionFC.style({
    color: '222222',
    fillColor: '00000000',
    width: 2
  }),
  {},
  '山西省边界',
  true
);

Map.addLayer(
  taiyuanPoint,
  {color: 'FF0000'},
  '太原核对点',
  false
);

print(
  '本步骤只请求一个年度的地图。' +
  '请等待地图加载，使用 Inspector 点击省内像元检查 GPP。' +
  '显示范围 0—2000 只是初始配色，不是数据有效值范围。'
);
//地图加载后，可以用右侧 Inspector 点击省内几个位置，确认 GPP 有数值；个别地点无值不一定是错误，可能是产品掩膜或本年完整性筛选造成的。

//这一步通过后，第二步再逐年导出年度 GPP 到 GEE Assets，保存中间结果；先不运行趋势统计。