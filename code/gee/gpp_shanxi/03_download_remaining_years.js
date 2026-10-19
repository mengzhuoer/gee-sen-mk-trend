/**************************************************************
 * 第三步：批量生成年度 GPP 下载链接（2001—2021）
 *
 * 数据：MODIS/061/MOD17A2HGF
 * 指标：全年各时段 GPP 累加，Int32 原始整数
 * 边界：GAUL 山西省（太原点位选取），与 2022 年一致
 * 分辨率：500 m（图像自身投影），与已下载的 2022 年方法一致
 * 比例：像元值 × 0.1 = g C/m²/year
 *
 * 2022 年已下载完成，本脚本生成 2001—2021 共 21 个链接。
 *
 * 用法：
 * 1. 点 Run，Console 按年份打印 21 条链接；
 * 2. 逐条点击下载（链接有有效期，建议一次点完）；
 * 3. 若某条过期（403），重新 Run 再点该年。
 *
 * 不做 Sen/MK、不生成图表、不创建导出任务。
 **************************************************************/


// ==================== 1. 参数 ====================

// 2022 已下载，这里生成其余年份。
var DOWNLOAD_START_YEAR = 2001;
var DOWNLOAD_END_YEAR = 2021;


// ==================== 2. 山西省边界 ====================

// 与第一步、第二步完全相同的选法：包含太原的省级行政区。
var taiyuanPoint = ee.Geometry.Point([112.549, 37.871]);

var regionFC = ee.FeatureCollection(
  'FAO/GAUL/2015/level1'
).filterBounds(taiyuanPoint);

var roi = regionFC.geometry();

print('第三步：批量生成年度 GPP 下载链接。');
print('年份范围：',
  DOWNLOAD_START_YEAR + '—' + DOWNLOAD_END_YEAR);
print('研究区要素数，预期为 1：', regionFC.size());
print('省级名称：', regionFC.aggregate_array('ADM1_NAME'));


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

print('本批次 8 天影像数量：', source.size());


// ==================== 4. 有效值处理 ====================

// 与第二步一致：产品原有掩膜保留，
// 有效范围 0—3000，真实 0 值保留，缺测不填 0。
function prepareGPP(image) {
  image = ee.Image(image);

  var raw = image.select('Gpp');

  return raw
    .updateMask(raw.gte(0).and(raw.lte(3000)))
    .rename('raw')
    .toInt32()
    .copyProperties(image, ['system:time_start']);
}


// ==================== 5. 单年累计（与第二步同一逻辑） ===

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

  // 本年必须有完整的 46 个时段。
  var yearIsComplete = ee.Image.constant(
    nImages.eq(46)
  ).setDefaultProjection(
    ee.Image(yearSource.first()).select('Gpp')
      .projection()
  );

  // 每个像元全部时段都有有效值才算完整年。
  var annualMask = validCount.eq(nImages)
    .and(yearIsComplete);

  var annualRaw = prepared.select('raw')
    .sum()
    .rename('GPP_raw_annual')
    .toInt32()
    .updateMask(annualMask)
    .clipToCollection(regionFC);

  return annualRaw.set({
    year: year,
    n_images: nImages,
    'system:time_start': start.millis()
  });
}


// ==================== 6. 逐年打印下载链接 ====================

// 客户端循环：每一年构建一次并打印链接。
// 构建链接本身不做全省统计，点链接时才触发该年计算。
for (
  var year = DOWNLOAD_START_YEAR;
  year <= DOWNLOAD_END_YEAR;
  year++
) {
  var annualRaw = buildAnnualRaw(year);

  var url = annualRaw.getDownloadURL({
    region: roi,
    scale: 500,
    format: 'GEO_TIFF',
    filename: 'Shanxi_GAUL_GPP_raw_annual_' + year
  });

  print(
    '===== ' + year + ' 年，下载后应为 2—3 MB ====='
  );
  print(url);
}

print(
  '共 ' +
  (DOWNLOAD_END_YEAR - DOWNLOAD_START_YEAR + 1) +
  ' 条链接已生成。逐条点击下载；' +
  '若某条过期，重新 Run 后只点对应年份即可。'
);
