/**************************************************************
 * 第五步：补缺 —— 下载 2017 和 2022（与已有文件同一网格）
 *
 * 背景：
 * - 2017：当初批量下载时漏点了一条链接；
 * - 2022：最早那版是 MODIS Sinusoidal 投影，与其余 20 个
 *   EPSG:4326 文件网格不一致，需要用同一口径重下。
 *
 * 本脚本一次 Run 打印两个年份：
 *   1) 该年的指纹（8 个点的取值，用于下载后核对身份）；
 *   2) 该年的下载链接。
 *
 * 点击两条链接分别下载（浏览器会存成哈希名），
 * 放进 D:\Zcode\GEE\导出 后，由分析端读指纹确认
 * 哪个文件是 2017、哪个是 2022，并改名归位。
 *
 * 网格/边界/口径与 03 脚本完全一致。
 **************************************************************/


// ==================== 1. 参数 ====================

// 需要补的两个年份
var FILL_YEARS = [2017, 2022];

// 与已有 20 个下载文件完全一致的网格
var GRID_CRS = 'EPSG:4326';
var GRID_TRANSFORM = [
  0.004491576420597608, 0, 110.22777693788589,
  0, -0.004491576420597608, 40.77902232260568
];

// 与 04 相同的 8 个指纹点
var POINTS = [
  [113.841250168, 39.613458241],
  [110.575874110, 34.749080978],
  [113.639129229, 36.222318044],
  [110.760028744, 38.095305411],
  [112.727339216, 38.104288564],
  [111.716734521, 36.235792773],
  [111.838007085, 39.586508783],
  [110.441126818, 36.743340909]
];


// ==================== 2. 山西省边界 ====================

var taiyuanPoint = ee.Geometry.Point([112.549, 37.871]);

var regionFC = ee.FeatureCollection(
  'FAO/GAUL/2015/level1'
).filterBounds(taiyuanPoint);

var roi = regionFC.geometry();

var pointsFC = ee.FeatureCollection(
  POINTS.map(function(p, idx) {
    return ee.Feature(ee.Geometry.Point(p), {
      point_index: idx
    });
  })
);


// ==================== 3. 原始数据 ====================

// 覆盖两个补缺年份的完整范围
var source = ee.ImageCollection(
  'MODIS/061/MOD17A2HGF'
)
  .filterDate(
    ee.Date.fromYMD(2017, 1, 1),
    ee.Date.fromYMD(2023, 1, 1)
  )
  .filterBounds(roi)
  .sort('system:time_start');

print('本批次影像数（2017 与 2022 共 92 张）：',
  source.size());


// ==================== 4. 单年累计（与下载脚本一致） =====

function prepareGPP(image) {
  image = ee.Image(image);

  var raw = image.select('Gpp');

  return raw
    .updateMask(raw.gte(0).and(raw.lte(3000)))
    .rename('raw')
    .toInt32()
    .copyProperties(image, ['system:time_start']);
}

function buildAnnualRaw(year) {
  var start = ee.Date.fromYMD(year, 1, 1);
  var end = start.advance(1, 'year');

  var yearSource = source.filterDate(start, end);
  var prepared = yearSource.map(prepareGPP);

  var nImages = yearSource.size();

  var validCount = prepared.select('raw')
    .count()
    .setDefaultProjection(
      ee.Image(yearSource.first()).select('Gpp')
        .projection()
    );

  var yearIsComplete = ee.Image.constant(
    nImages.eq(46)
  ).setDefaultProjection(
    ee.Image(yearSource.first()).select('Gpp')
      .projection()
  );

  var annualMask = validCount.eq(nImages)
    .and(yearIsComplete);

  return prepared.select('raw')
    .sum()
    .rename('GPP_raw_annual')
    .toInt32()
    .updateMask(annualMask)
    .clipToCollection(regionFC)
    .set({year: year, n_images: nImages});
}


// ==================== 5. 打印指纹 + 下载链接 ====================

FILL_YEARS.forEach(function(year) {
  year = ee.Number(year);

  var annual = buildAnnualRaw(year);

  // 指纹：与 04 同一网格同一批点，用于下载后核对文件身份
  var snap = annual.reproject(GRID_CRS, GRID_TRANSFORM);
  var res = snap.reduceRegions({
    collection: pointsFC,
    reducer: ee.Reducer.first()
  });
  var vals = res.aggregate_array('first').getInfo();

  print('===== ' + year.getInfo() + ' 年 =====');
  print('指纹（下载后核对用）: [' + String(vals) + ']');

  // 时段数量检查，预期 46
  print('时段数量（应为 46）:', annual.get('n_images'));

  // 下载链接
  var url = annual.getDownloadURL({
    region: roi,
    scale: 500,
    format: 'GEO_TIFF'
  });
  print(url);
});

print(
  '两条链接分别对应上面打印指纹的年份（按出现顺序）。' +
  '下载后把文件放进 D:\\Zcode\\GEE\\导出，' +
  '由分析端核对指纹并改名。链接过期就重新 Run。'
);
