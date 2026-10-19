/**************************************************************
 * 第四步：年度指纹 —— 确定每个下载文件对应哪一年
 *
 * 背景：
 * getDownloadURL 返回的文件名是哈希名，文件内没有年份标签，
 * 无法直接知道每个 tiff 是哪一年。解决办法：
 *
 * 1. 本地在 20 个文件的同一组像元中心取 8 个点，得到 20 组指纹；
 * 2. 本脚本在与下载文件完全一致的网格上（EPSG:4326 + 同一
 *    geotransform）对 2001—2021 每年在同样的 8 个点取值并打印；
 * 3. 两边的 8 元素向量一一比对，即可把文件名映射到年份。
 *
 * 网格参数直接取自已下载文件的 geotransform，reproject 保证
 * 采样网格与文件逐像元一致。
 *
 * 输出格式（逐行，复制全部行粘贴给分析端即可）：
 *   2001 => [5880, 4798, ...]
 *   ...
 *
 * 不做趋势分析、不生成图表、不下载文件。
 **************************************************************/


// ==================== 1. 参数 ====================

var DOWNLOAD_START_YEAR = 2001;
var DOWNLOAD_END_YEAR = 2021;

// 与已下载 4326 文件完全一致的网格
// （从本地 tiff 的 geotransform 原样读出）
var GRID_CRS = 'EPSG:4326';
var GRID_TRANSFORM = [
  0.004491576420597608, 0, 110.22777693788589,
  0, -0.004491576420597608, 40.77902232260568
];

// 8 个采样点：下载网格的像元中心，省内，
// 已验证在 20 个本地文件上全部为有效值。
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

var source = ee.ImageCollection(
  'MODIS/061/MOD17A2HGF'
)
  .filterDate(
    ee.Date.fromYMD(DOWNLOAD_START_YEAR, 1, 1),
    ee.Date.fromYMD(DOWNLOAD_END_YEAR + 1, 1, 1)
  )
  .filterBounds(roi)
  .sort('system:time_start');

print('指纹计算开始。影像数（应为 21 年 × 46 = 966）：',
  source.size());


// ==================== 4. 单年累计（与下载脚本完全一致） =

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
    .clipToCollection(regionFC);
}


// ==================== 5. 逐年打印指纹 ====================

print('--- FINGERPRINT START ---');

for (
  var year = DOWNLOAD_START_YEAR;
  year <= DOWNLOAD_END_YEAR;
  year++
) {
  var annual = buildAnnualRaw(year);

  // 重投影到下载文件的精确网格，保证逐像元一致
  var snap = annual.reproject(GRID_CRS, GRID_TRANSFORM);

  var res = snap.reduceRegions({
    collection: pointsFC,
    reducer: ee.Reducer.first()
  });

  // 客户端取值，打印成单行，方便整体复制
  var vals = res.aggregate_array('first').getInfo();

  print(year + ' => [' + String(vals) + ']');
}

print('--- FINGERPRINT END ---');
print(
  '请把 START 到 END 之间的所有行复制粘贴。' +
  '本步骤只做点采样，不生成任何图表。'
);
