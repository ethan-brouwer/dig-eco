/*
  FILE: gee/turbidity/san_sebastian_impact_pool_image_turbidity.js
  PURPOSE: Build a per-image Sentinel-2 turbidity time series for the impact_pool polygon.

  GEE IMPORTS REQUIRED
  - impact_pool
*/

// ================= USER SETTINGS =================
var siteName = "San Sebastian impact pool";
var analysisStartYear = 2017;
var cloudMax = 80;
var analysisScaleMeters = 10;
var useWetSeasonFilter = false;
var wetSeasonMonths = [5, 6, 7, 8, 9, 10];
var ndwiThreshold = -0.05;
var mndwiThreshold = -0.05;
var nirWaterMax = 0.15;
var swirWaterMax = 0.12;
var minWaterPixels = 3;
var noDataValue = -9999;

// ================= GEOMETRY HELPERS =================
function toGeometry(obj) {
  var t = ee.String(ee.Algorithms.ObjectType(obj));
  return ee.Geometry(ee.Algorithms.If(
    t.compareTo("FeatureCollection").eq(0),
    ee.FeatureCollection(obj).geometry(),
    ee.Algorithms.If(
      t.compareTo("Feature").eq(0),
      ee.Feature(obj).geometry(),
      ee.Geometry(obj)
    )
  ));
}

function safeNumber(x) {
  return ee.Number(ee.Algorithms.If(ee.Algorithms.IsEqual(x, null), noDataValue, x));
}

var impactPoolGeom = toGeometry(impact_pool);
var aoi = impactPoolGeom.bounds().buffer(200);

// ================= MAP SETUP =================
Map.setOptions("SATELLITE");
Map.centerObject(impactPoolGeom, 16);
Map.addLayer(impactPoolGeom, {color: "FFFF00"}, "impact_pool", true);
Map.addLayer(aoi, {color: "FFFFFF"}, "AOI", false);

// ================= DATE WINDOW =================
var endDate = ee.Date(Date.now());
var startDate = ee.Date.fromYMD(analysisStartYear, 1, 1);

print("Analysis window", startDate.format("YYYY-MM-dd"), endDate.format("YYYY-MM-dd"));
print("Scale (m)", analysisScaleMeters);
print("Wet-season month filter", useWetSeasonFilter ? wetSeasonMonths : "OFF");
print("NDWI threshold", ndwiThreshold);
print("MNDWI threshold", mndwiThreshold);
print("NIR max", nirWaterMax);
print("SWIR1 max", swirWaterMax);
print("Min water pixels", minWaterPixels);

// ================= S2 HELPERS =================
function maskS2(img) {
  var scl = img.select("SCL");
  var qa = img.select("QA60");

  var cloud = qa.bitwiseAnd(1 << 10).neq(0)
    .or(qa.bitwiseAnd(1 << 11).neq(0));
  var shadow = scl.eq(3);
  var cirrus = scl.eq(10);
  var snow = scl.eq(11);
  var saturated = scl.eq(1);

  var bad = cloud.or(shadow).or(cirrus).or(snow).or(saturated);
  return img.updateMask(bad.not());
}

function scaleS2(img) {
  var optical = img.select(["B2", "B3", "B4", "B8", "B11"]).multiply(0.0001);
  return img.addBands(optical, null, true);
}

function addMonth(img) {
  return img.set("month", ee.Date(img.get("system:time_start")).get("month"));
}

function addIndices(img) {
  var ndwi = img.normalizedDifference(["B3", "B8"]).rename("NDWI");
  var mndwi = img.normalizedDifference(["B3", "B11"]).rename("MNDWI");
  var nsmi = img.expression(
    "(RED + GREEN - BLUE) / (RED + GREEN + BLUE)",
    {
      RED: img.select("B4"),
      GREEN: img.select("B3"),
      BLUE: img.select("B2")
    }
  ).rename("NSMI");
  return img.addBands([ndwi, mndwi, nsmi]);
}

function waterMask(img) {
  return img.select("NDWI").gt(ndwiThreshold)
    .or(img.select("MNDWI").gt(mndwiThreshold))
    .or(
      img.select("B8").lt(nirWaterMax)
        .and(img.select("B11").lt(swirWaterMax))
    )
    .selfMask()
    .rename("WATER_MASK");
}

var s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  .filterBounds(aoi)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", cloudMax))
  .map(maskS2)
  .map(scaleS2)
  .map(addMonth)
  .map(addIndices);

if (useWetSeasonFilter) {
  s2 = s2.filter(ee.Filter.inList("month", wetSeasonMonths));
}

print("Raw image count", s2.size());

// ================= PER-IMAGE TIME SERIES =================
var scoredS2 = s2.map(function(img) {
  var date = ee.Date(img.get("system:time_start"));
  var water = waterMask(img).clip(impactPoolGeom);
  var validObs = img.select("B4").mask().clip(impactPoolGeom);

  var validPx = safeNumber(validObs.reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: impactPoolGeom,
    scale: analysisScaleMeters,
    bestEffort: true,
    maxPixels: 1e9
  }).get("B4"));

  var waterPx = safeNumber(water.reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: impactPoolGeom,
    scale: analysisScaleMeters,
    bestEffort: true,
    maxPixels: 1e9
  }).get("WATER_MASK"));

  var waterOnly = img.select(["NDWI", "MNDWI", "NSMI"]).updateMask(water);
  var stats = waterOnly.reduceRegion({
    reducer: ee.Reducer.mean()
      .combine({reducer2: ee.Reducer.median(), sharedInputs: true})
      .combine({reducer2: ee.Reducer.stdDev(), sharedInputs: true}),
    geometry: impactPoolGeom,
    scale: analysisScaleMeters,
    bestEffort: true,
    maxPixels: 1e9
  });

  var isValid = waterPx.gte(minWaterPixels);

  return img.set({
    image_id: ee.String(img.get("PRODUCT_ID")),
    granule_id: ee.String(img.get("GRANULE_ID")),
    date: date.format("YYYY-MM-dd"),
    year: ee.Number.parse(date.format("YYYY")),
    month_num: ee.Number.parse(date.format("M")),
    cloudy_pixel_pct: img.get("CLOUDY_PIXEL_PERCENTAGE"),
    valid_px: validPx,
    water_px: waterPx,
    qa_flag: ee.Algorithms.If(isValid, "water_mask", "low_water_px"),
    ndwi_mean: ee.Algorithms.If(isValid, safeNumber(stats.get("NDWI_mean")), noDataValue),
    ndwi_median: ee.Algorithms.If(isValid, safeNumber(stats.get("NDWI_median")), noDataValue),
    ndwi_stddev: ee.Algorithms.If(isValid, safeNumber(stats.get("NDWI_stdDev")), noDataValue),
    mndwi_mean: ee.Algorithms.If(isValid, safeNumber(stats.get("MNDWI_mean")), noDataValue),
    mndwi_median: ee.Algorithms.If(isValid, safeNumber(stats.get("MNDWI_median")), noDataValue),
    mndwi_stddev: ee.Algorithms.If(isValid, safeNumber(stats.get("MNDWI_stdDev")), noDataValue),
    nsmi_mean: ee.Algorithms.If(isValid, safeNumber(stats.get("NSMI_mean")), noDataValue),
    nsmi_median: ee.Algorithms.If(isValid, safeNumber(stats.get("NSMI_median")), noDataValue),
    nsmi_stddev: ee.Algorithms.If(isValid, safeNumber(stats.get("NSMI_stdDev")), noDataValue)
  });
});

var validScoredS2 = scoredS2.filter(ee.Filter.eq("qa_flag", "water_mask"));

var imageStats = ee.FeatureCollection(scoredS2.map(function(img) {
  return ee.Feature(null, img.toDictionary([
    "system:time_start",
    "image_id",
    "granule_id",
    "date",
    "year",
    "month_num",
    "cloudy_pixel_pct",
    "valid_px",
    "water_px",
    "qa_flag",
    "ndwi_mean",
    "ndwi_median",
    "ndwi_stddev",
    "mndwi_mean",
    "mndwi_median",
    "mndwi_stddev",
    "nsmi_mean",
    "nsmi_median",
    "nsmi_stddev"
  ]));
})).sort("system:time_start");

var validImageStats = imageStats.filter(ee.Filter.eq("qa_flag", "water_mask"));

print("Per-image turbidity table", imageStats);
print("Valid image count", validImageStats.size());

var highNsmiEvents = validImageStats.sort("nsmi_mean", false).limit(15);
var lowNdwiEvents = validImageStats.sort("ndwi_mean", true).limit(15);
var lowMndwiEvents = validImageStats.sort("mndwi_mean", true).limit(15);

print("Highest NSMI images", highNsmiEvents);
print("Lowest NDWI images", lowNdwiEvents);
print("Lowest MNDWI images", lowMndwiEvents);

// ================= CHARTS =================
var chartFc = imageStats.map(function(f) {
  function cleanValue(name) {
    var value = f.get(name);
    return ee.Algorithms.If(ee.Algorithms.IsEqual(value, noDataValue), null, value);
  }

  return f.set({
    ndwi_mean_chart: cleanValue("ndwi_mean"),
    mndwi_mean_chart: cleanValue("mndwi_mean"),
    nsmi_mean_chart: cleanValue("nsmi_mean"),
    water_px_chart: ee.Algorithms.If(ee.Number(f.get("water_px")).gt(0), f.get("water_px"), null)
  });
});

var ndwiChart = ui.Chart.feature.byFeature(
  chartFc.filter(ee.Filter.notNull(["ndwi_mean_chart"])),
  "date",
  ["ndwi_mean_chart", "mndwi_mean_chart"]
).setChartType("LineChart")
  .setOptions({
    title: "Impact Pool Water Signal: NDWI and MNDWI",
    hAxis: {title: "Image date", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "Index value"},
    lineWidth: 2,
    pointSize: 3,
    colors: ["#1E88E5", "#43A047"]
  });
print(ndwiChart);

var nsmiChart = ui.Chart.feature.byFeature(
  chartFc.filter(ee.Filter.notNull(["nsmi_mean_chart"])),
  "date",
  ["nsmi_mean_chart"]
).setChartType("LineChart")
  .setOptions({
    title: "Impact Pool Water Signal: NSMI",
    hAxis: {title: "Image date", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "NSMI mean"},
    lineWidth: 2,
    pointSize: 3,
    colors: ["#E53935"]
  });
print(nsmiChart);

var waterPxChart = ui.Chart.feature.byFeature(
  chartFc.filter(ee.Filter.notNull(["water_px_chart"])),
  "date",
  ["water_px_chart"]
).setChartType("ColumnChart")
  .setOptions({
    title: "Impact Pool Water Pixels Used Per Image",
    hAxis: {title: "Image date", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "Water pixels at 10 m"},
    legend: {position: "none"},
    colors: ["#26A69A"]
  });
print(waterPxChart);

// ================= QUICKLOOKS =================
var latestValid = ee.Image(ee.Algorithms.If(
  validScoredS2.size().gt(0),
  validScoredS2.sort("system:time_start", false).first(),
  s2.sort("system:time_start", false).first()
));
print("Preview image source", ee.Algorithms.If(
  validScoredS2.size().gt(0),
  "latest valid image",
  "latest raw image"
));
Map.addLayer(
  latestValid.select(["B4", "B3", "B2"]).clip(impactPoolGeom),
  {min: 0.02, max: 0.3},
  "Latest valid RGB",
  false
);
Map.addLayer(
  latestValid.select("NSMI").clip(impactPoolGeom),
  {min: -0.3, max: 0.3, palette: ["#2166AC", "#F7F7F7", "#B2182B"]},
  "Latest valid NSMI",
  false
);
