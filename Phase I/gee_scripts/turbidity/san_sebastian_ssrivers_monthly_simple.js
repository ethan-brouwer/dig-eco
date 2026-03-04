/*
  FILE: Phase I/gee_scripts/turbidity/san_sebastian_ssrivers_monthly_simple.js
  PURPOSE: Simple SSrivers-buffer workflow for water proxy mapping + monthly charts.

  REQUIREMENT
  - In GEE Imports, provide `SSrivers` as a Geometry, Feature, or FeatureCollection.
*/

// ================= USER SETTINGS =================
var siteName = "San Sebastian Mine (MRDS)";
var siteLon = -87.93002;
var siteLat = 13.6509;

var recentMonths = 12;
var cloudMax = 60;
var corridorBufferMeters = 120;
var useWetSeasonFilter = false;
var wetSeasonMonths = [5, 6, 7, 8, 9, 10];

// Water mask thresholds (keep simple and easy to tune)
var ndwiThreshold = -0.02;
var mndwiThreshold = -0.05;
var minObsCount = 2;
var noDataValue = -9999;

// ================= MAP SETUP =================
Map.setOptions("SATELLITE");
var sitePoint = ee.Geometry.Point([siteLon, siteLat]);
Map.centerObject(sitePoint, 12);
Map.addLayer(sitePoint, {color: "FF0000"}, siteName, true);

// ================= IMPORTED GEOMETRY (SSrivers) =================
var ssType = ee.String(ee.Algorithms.ObjectType(SSrivers));
var ssGeom = ee.Geometry(ee.Algorithms.If(
  ssType.compareTo("FeatureCollection").eq(0),
  ee.FeatureCollection(SSrivers).geometry(),
  ee.Algorithms.If(
    ssType.compareTo("Feature").eq(0),
    ee.Feature(SSrivers).geometry(),
    ee.Geometry(SSrivers)
  )
));

var riverCorridor = ssGeom.buffer(corridorBufferMeters);
var aoi = riverCorridor.bounds().buffer(1000);
var corridorMask = ee.Image.constant(1).clip(riverCorridor).selfMask().rename("CORRIDOR_MASK");

Map.addLayer(ssGeom, {color: "00E5FF"}, "SSrivers line", true);
Map.addLayer(riverCorridor, {color: "FFD54F"}, "SSrivers corridor", true);
Map.addLayer(aoi, {color: "FFFFFF"}, "AOI", false);

// ================= DATE WINDOW =================
var endDate = ee.Date(Date.now());
var startDate = endDate.advance(-recentMonths, "month");
print("Recent window", startDate.format("YYYY-MM-dd"), endDate.format("YYYY-MM-dd"));
print("Corridor buffer (m)", corridorBufferMeters);
print("Wet season filter", useWetSeasonFilter ? wetSeasonMonths : "OFF");
print("NDWI threshold", ndwiThreshold);
print("MNDWI threshold", mndwiThreshold);

// ================= S2 HELPERS =================
function maskS2(img) {
  var scl = img.select("SCL");
  var qa = img.select("QA60");
  var cloud = qa.bitwiseAnd(1 << 10).neq(0).or(qa.bitwiseAnd(1 << 11).neq(0));
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
  var green = img.select("B3");
  var red = img.select("B4");
  var ndwi = img.normalizedDifference(["B3", "B8"]).rename("NDWI");
  var mndwi = img.normalizedDifference(["B3", "B11"]).rename("MNDWI");
  var ndti = img.normalizedDifference(["B4", "B3"]).rename("NDTI");
  var tssProxy = red.multiply(1000).rename("TSS_PROXY");
  var redGreen = red.divide(green).rename("RED_GREEN");
  return img.addBands([ndwi, mndwi, ndti, tssProxy, redGreen]);
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

print("Image count", s2.size());

// ================= CURRENT MAP LAYERS =================
var recent = s2.median().clip(aoi);
var obsCount = s2.select("B4").count().clip(aoi).rename("OBS_COUNT");
var enoughObs = obsCount.gte(minObsCount);

var riverMask = recent.select("NDWI").gt(ndwiThreshold)
  .or(recent.select("MNDWI").gt(mndwiThreshold))
  .updateMask(enoughObs)
  .updateMask(corridorMask)
  .selfMask()
  .rename("RIVER_MASK");

var tssWater = recent.select("TSS_PROXY").updateMask(riverMask).rename("TSS_WATER");
var ndtiWater = recent.select("NDTI").updateMask(riverMask).rename("NDTI_WATER");
var redGreenWater = recent.select("RED_GREEN").updateMask(riverMask).rename("RED_GREEN_WATER");

var waterProxyComposite = ee.Image.cat([
  recent.select("TSS_PROXY").unitScale(20, 220).clamp(0, 1),
  recent.select("NDTI").unitScale(-0.1, 0.25).clamp(0, 1),
  recent.select("RED_GREEN").unitScale(0.7, 2.0).clamp(0, 1)
]).updateMask(riverMask).rename(["TSS_RGB", "NDTI_RGB", "RED_GREEN_RGB"]);

Map.addLayer(recent, {bands: ["B4", "B3", "B2"], min: 0.02, max: 0.30}, "Recent composite (true color)", true);
Map.addLayer(recent, {bands: ["B8", "B4", "B3"], min: 0.03, max: 0.45}, "Recent composite (false color)", false);
Map.addLayer(obsCount, {min: 0, max: 30, palette: ["2b2b2b", "f7f7f7", "00ff00"]}, "QA observation count", false);
Map.addLayer(corridorMask, {palette: ["FFF59D"]}, "Corridor mask", false);
Map.addLayer(riverMask, {palette: ["FFD54F"]}, "River mask", true);
Map.addLayer(waterProxyComposite, {bands: ["TSS_RGB", "NDTI_RGB", "RED_GREEN_RGB"], min: 0, max: 1}, "Water proxy composite", true);
Map.addLayer(tssWater, {min: 20, max: 300, palette: ["08306b", "41b6c4", "ffffbf", "fdae61", "d73027"]}, "TSS proxy (water)", false);
Map.addLayer(ndtiWater, {min: -0.3, max: 0.4, palette: ["2166AC", "F7F7F7", "B2182B"]}, "NDTI (water)", false);
Map.addLayer(redGreenWater, {min: 0.7, max: 2.0, palette: ["F7FBFF", "FDAE61", "D73027"]}, "Red/Green ratio (water)", false);

// ================= MONTHLY GRAPH OUTPUT =================
function monthStartList(start, end) {
  var monthCount = ee.Number(end.difference(start, "month")).floor();
  return ee.List.sequence(0, monthCount.subtract(1)).map(function(m) {
    return start.advance(ee.Number(m), "month");
  });
}

function safeNumber(x) {
  return ee.Number(ee.Algorithms.If(ee.Algorithms.IsEqual(x, null), noDataValue, x));
}

var monthlyStarts = monthStartList(
  ee.Date(startDate.format("YYYY-MM-01")),
  ee.Date(endDate.format("YYYY-MM-01")).advance(1, "month")
);

var monthlyStats = ee.FeatureCollection(monthlyStarts.map(function(mStart) {
  mStart = ee.Date(mStart);
  var mEnd = mStart.advance(1, "month");
  var col = s2.filterDate(mStart, mEnd);
  var count = col.size();

  var emptyFeature = ee.Feature(null, {
    "system:time_start": mStart.millis(),
    date: mStart.format("YYYY-MM"),
    image_count: count,
    water_px: 0,
    tss_proxy: noDataValue,
    ndti: noDataValue,
    red_green: noDataValue
  });

  return ee.Feature(ee.Algorithms.If(count.gt(0), (function() {
    var img = ee.Image(col.median()).clip(aoi);
    var mObs = col.select("B4").count().clip(aoi);
    var mMask = img.select("NDWI").gt(ndwiThreshold)
      .or(img.select("MNDWI").gt(mndwiThreshold))
      .updateMask(mObs.gte(minObsCount))
      .updateMask(corridorMask)
      .selfMask();

    var waterPx = safeNumber(mMask.reduceRegion({
      reducer: ee.Reducer.count(),
      geometry: riverCorridor,
      scale: 20,
      bestEffort: true,
      maxPixels: 1e9
    }).values().get(0));

    var stats = img.select(["TSS_PROXY", "NDTI", "RED_GREEN"])
      .updateMask(mMask)
      .reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: riverCorridor,
        scale: 20,
        bestEffort: true,
        maxPixels: 1e9
      });

    var enoughWater = waterPx.gte(5);

    return ee.Feature(null, {
      "system:time_start": mStart.millis(),
      date: mStart.format("YYYY-MM"),
      image_count: count,
      water_px: waterPx,
      tss_proxy: ee.Algorithms.If(enoughWater, safeNumber(stats.get("TSS_PROXY")), noDataValue),
      ndti: ee.Algorithms.If(enoughWater, safeNumber(stats.get("NDTI")), noDataValue),
      red_green: ee.Algorithms.If(enoughWater, safeNumber(stats.get("RED_GREEN")), noDataValue)
    });
  })(), emptyFeature));
})).sort("system:time_start");

print("Monthly table", monthlyStats);

var tssChart = ui.Chart.feature.byFeature(monthlyStats, "date", ["tss_proxy"])
  .setChartType("LineChart")
  .setOptions({
    title: "Monthly TSS Proxy (SSrivers Corridor)",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "TSS proxy"},
    lineWidth: 2,
    pointSize: 3,
    colors: ["#D84315"]
  });
print(tssChart);

var ndtiChart = ui.Chart.feature.byFeature(monthlyStats, "date", ["ndti"])
  .setChartType("LineChart")
  .setOptions({
    title: "Monthly NDTI (SSrivers Corridor)",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "NDTI"},
    lineWidth: 2,
    pointSize: 3,
    colors: ["#8E24AA"]
  });
print(ndtiChart);

var redGreenChart = ui.Chart.feature.byFeature(monthlyStats, "date", ["red_green"])
  .setChartType("LineChart")
  .setOptions({
    title: "Monthly Red/Green Ratio (SSrivers Corridor)",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "Red/Green"},
    lineWidth: 2,
    pointSize: 3,
    colors: ["#2E7D32"]
  });
print(redGreenChart);
