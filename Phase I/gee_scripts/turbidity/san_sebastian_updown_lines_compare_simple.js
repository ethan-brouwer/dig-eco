/*
  FILE: Phase I/gee_scripts/turbidity/san_sebastian_updown_lines_compare_simple.js
  PURPOSE: Compare upstream1 vs downstream1 proxy behavior by month.

  REQUIREMENT (GEE Imports)
  - impact_point
  - upstream1
  - downstream1
*/

// ================= USER SETTINGS =================
var siteName = "San Sebastian Mine (MRDS)";
var siteLon = -87.93002;
var siteLat = 13.6509;

var analysisStartYear = 2020;
var cloudMax = 60;
var corridorBufferMeters = 50;
var useWetSeasonFilter = true;
var wetSeasonMonths = [5, 6, 7, 8, 9, 10];

var ndwiThreshold = -0.02;
var mndwiThreshold = -0.05;
var minObsCount = 2;
var minAreaPixelsPerMonth = 5;
var useCorridorDirectly = true;
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

var impactGeom = toGeometry(impact_point);
var upstreamGeom = toGeometry(upstream1);
var downstreamGeom = toGeometry(downstream1);

var upstreamCorridor = upstreamGeom.buffer(corridorBufferMeters);
var downstreamCorridor = downstreamGeom.buffer(corridorBufferMeters);
var combinedCorridor = upstreamCorridor.union(downstreamCorridor, 1);
var aoi = combinedCorridor.bounds().buffer(1000);

var reaches = ee.FeatureCollection([
  ee.Feature(upstreamCorridor, {reach_id: "upstream1", reach_type: "upstream"}),
  ee.Feature(downstreamCorridor, {reach_id: "downstream1", reach_type: "downstream"})
]);

// ================= MAP SETUP =================
Map.setOptions("SATELLITE");
Map.centerObject(impactGeom, 12);
Map.addLayer(ee.Geometry.Point([siteLon, siteLat]), {color: "FF0000"}, siteName, false);
Map.addLayer(impactGeom, {color: "FFFF00"}, "impact_point", true);
Map.addLayer(upstreamGeom, {color: "2962FF"}, "upstream1 line", true);
Map.addLayer(downstreamGeom, {color: "FF6D00"}, "downstream1 line", true);
Map.addLayer(upstreamCorridor, {color: "2962FF"}, "upstream corridor", false);
Map.addLayer(downstreamCorridor, {color: "FF6D00"}, "downstream corridor", false);
Map.addLayer(aoi, {color: "FFFFFF"}, "AOI", false);

// ================= DATE WINDOW =================
var currentYear = new Date().getUTCFullYear();
var startYear = analysisStartYear;
var endDate = ee.Date(Date.now());
var startDate = ee.Date.fromYMD(startYear, 1, 1);

print("Analysis window", startDate.format("YYYY-MM-dd"), endDate.format("YYYY-MM-dd"));
print("Analysis start year", analysisStartYear);
print("Buffer (m)", corridorBufferMeters);
print("Wet season filter", useWetSeasonFilter ? wetSeasonMonths : "OFF");
print("NDWI threshold", ndwiThreshold);
print("MNDWI threshold", mndwiThreshold);
print("Use corridor directly", useCorridorDirectly);

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

function corridorMask(geom) {
  return ee.Image.constant(1).clip(geom).selfMask();
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

// ================= CURRENT VISUALIZATION =================
var recent = s2.median().clip(aoi);
var obsCount = s2.select("B4").count().clip(aoi);
var enoughObs = obsCount.gte(minObsCount);

var upMaskCurrent = recent.select("NDWI").gt(ndwiThreshold)
  .or(recent.select("MNDWI").gt(mndwiThreshold))
  .updateMask(enoughObs)
  .updateMask(corridorMask(upstreamCorridor))
  .selfMask();

var downMaskCurrent = recent.select("NDWI").gt(ndwiThreshold)
  .or(recent.select("MNDWI").gt(mndwiThreshold))
  .updateMask(enoughObs)
  .updateMask(corridorMask(downstreamCorridor))
  .selfMask();

var currentMask = upMaskCurrent.or(downMaskCurrent).selfMask().rename("CURRENT_MASK");
var waterProxyComposite = ee.Image.cat([
  recent.select("TSS_PROXY").unitScale(20, 220).clamp(0, 1),
  recent.select("NDTI").unitScale(-0.1, 0.25).clamp(0, 1),
  recent.select("RED_GREEN").unitScale(0.7, 2.0).clamp(0, 1)
]).updateMask(currentMask).rename(["TSS_RGB", "NDTI_RGB", "RED_GREEN_RGB"]);

Map.addLayer(recent, {bands: ["B4", "B3", "B2"], min: 0.02, max: 0.30}, "Recent composite (true color)", true);
Map.addLayer(recent, {bands: ["B8", "B4", "B3"], min: 0.03, max: 0.45}, "Recent composite (false color)", false);
Map.addLayer(obsCount, {min: 0, max: 30, palette: ["2b2b2b", "f7f7f7", "00ff00"]}, "QA observation count", false);
Map.addLayer(upMaskCurrent, {palette: ["2962FF"]}, "Current upstream mask", false);
Map.addLayer(downMaskCurrent, {palette: ["FF6D00"]}, "Current downstream mask", false);
Map.addLayer(currentMask, {palette: ["FFD54F"]}, "Current combined mask", true);
Map.addLayer(waterProxyComposite, {bands: ["TSS_RGB", "NDTI_RGB", "RED_GREEN_RGB"], min: 0, max: 1}, "Water proxy composite", true);

// ================= MONTHLY STATS =================
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

var monthlyLong = ee.FeatureCollection(monthlyStarts.map(function(mStart) {
  mStart = ee.Date(mStart);
  var mEnd = mStart.advance(1, "month");
  var col = s2.filterDate(mStart, mEnd);
  var count = col.size();

  return reaches.map(function(reach) {
    reach = ee.Feature(reach);
    var rGeom = reach.geometry();
    var rMask = corridorMask(rGeom);

    var baseFeature = ee.Feature(null, {
      "system:time_start": mStart.millis(),
      date: mStart.format("YYYY-MM"),
      year: ee.Number.parse(mStart.format("YYYY")),
      month_num: ee.Number.parse(mStart.format("M")),
      month_label: mStart.format("MMM"),
      reach_id: reach.get("reach_id"),
      reach_type: reach.get("reach_type"),
      image_count: count,
      water_px: 0,
      qa_flag: "no_images",
      tss_proxy: noDataValue,
      ndti: noDataValue,
      red_green: noDataValue
    });

    return ee.Feature(ee.Algorithms.If(count.gt(0), (function() {
      var img = ee.Image(col.median()).clip(aoi);
      var mObs = col.select("B4").count().clip(aoi);
      var monthWaterMask = img.select("NDWI").gt(ndwiThreshold)
        .or(img.select("MNDWI").gt(mndwiThreshold))
        .updateMask(mObs.gte(minObsCount))
        .updateMask(rMask)
        .selfMask();

      var monthAreaMask = rMask.updateMask(mObs.gte(minObsCount)).selfMask();
      var analysisMask = ee.Image(ee.Algorithms.If(useCorridorDirectly, monthAreaMask, monthWaterMask));

      var waterPx = safeNumber(analysisMask.reduceRegion({
        reducer: ee.Reducer.count(),
        geometry: rGeom,
        scale: 20,
        bestEffort: true,
        maxPixels: 1e9
      }).values().get(0));

      var stats = img.select(["TSS_PROXY", "NDTI", "RED_GREEN"])
        .updateMask(analysisMask)
        .reduceRegion({
          reducer: ee.Reducer.mean(),
          geometry: rGeom,
          scale: 20,
          bestEffort: true,
          maxPixels: 1e9
        });

      var enoughArea = waterPx.gte(minAreaPixelsPerMonth);
      return baseFeature.set({
        water_px: waterPx,
        qa_flag: ee.Algorithms.If(useCorridorDirectly, "corridor_direct", "water_mask"),
        tss_proxy: ee.Algorithms.If(enoughArea, safeNumber(stats.get("TSS_PROXY")), noDataValue),
        ndti: ee.Algorithms.If(enoughArea, safeNumber(stats.get("NDTI")), noDataValue),
        red_green: ee.Algorithms.If(enoughArea, safeNumber(stats.get("RED_GREEN")), noDataValue)
      });
    })(), baseFeature));
  });
})).flatten().sort("system:time_start");

print("Monthly long table", monthlyLong);

var monthlyChartFc = monthlyLong.map(function(f) {
  var tss = f.get("tss_proxy");
  var ndti = f.get("ndti");
  var rg = f.get("red_green");
  return f.set({
    tss_chart: ee.Algorithms.If(ee.Algorithms.IsEqual(tss, noDataValue), null, tss),
    ndti_chart: ee.Algorithms.If(ee.Algorithms.IsEqual(ndti, noDataValue), null, ndti),
    red_green_chart: ee.Algorithms.If(ee.Algorithms.IsEqual(rg, noDataValue), null, rg)
  });
});

var tssCompareChart = ui.Chart.feature.groups(monthlyChartFc.filter(ee.Filter.notNull(["tss_chart"])), "date", "tss_chart", "reach_type")
  .setChartType("LineChart")
  .setOptions({
    title: "Monthly TSS Proxy: Upstream vs Downstream",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "TSS proxy"},
    lineWidth: 2,
    pointSize: 3,
    colors: ["#2962FF", "#FF6D00"]
  });
print(tssCompareChart);

var ndtiCompareChart = ui.Chart.feature.groups(monthlyChartFc.filter(ee.Filter.notNull(["ndti_chart"])), "date", "ndti_chart", "reach_type")
  .setChartType("LineChart")
  .setOptions({
    title: "Monthly NDTI: Upstream vs Downstream",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "NDTI"},
    lineWidth: 2,
    pointSize: 3,
    colors: ["#2962FF", "#FF6D00"]
  });
print(ndtiCompareChart);

var redGreenCompareChart = ui.Chart.feature.groups(monthlyChartFc.filter(ee.Filter.notNull(["red_green_chart"])), "date", "red_green_chart", "reach_type")
  .setChartType("LineChart")
  .setOptions({
    title: "Monthly Red/Green Ratio: Upstream vs Downstream",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "Red/Green"},
    lineWidth: 2,
    pointSize: 3,
    colors: ["#2962FF", "#FF6D00"]
  });
print(redGreenCompareChart);

var monthlyDiff = ee.FeatureCollection(monthlyStarts.map(function(mStart) {
  mStart = ee.Date(mStart);
  var key = mStart.format("YYYY-MM");
  var monthFc = monthlyLong.filter(ee.Filter.eq("date", key));
  var up = ee.Feature(monthFc.filter(ee.Filter.eq("reach_type", "upstream")).first());
  var down = ee.Feature(monthFc.filter(ee.Filter.eq("reach_type", "downstream")).first());

  var upTss = ee.Number(up.get("tss_proxy"));
  var downTss = ee.Number(down.get("tss_proxy"));
  var upNdti = ee.Number(up.get("ndti"));
  var downNdti = ee.Number(down.get("ndti"));
  var upRg = ee.Number(up.get("red_green"));
  var downRg = ee.Number(down.get("red_green"));

  var tssDiff = ee.Algorithms.If(upTss.eq(noDataValue).or(downTss.eq(noDataValue)), noDataValue, downTss.subtract(upTss));
  var ndtiDiff = ee.Algorithms.If(upNdti.eq(noDataValue).or(downNdti.eq(noDataValue)), noDataValue, downNdti.subtract(upNdti));
  var rgDiff = ee.Algorithms.If(upRg.eq(noDataValue).or(downRg.eq(noDataValue)), noDataValue, downRg.subtract(upRg));

  return ee.Feature(null, {
    "system:time_start": mStart.millis(),
    date: key,
    tss_diff: tssDiff,
    ndti_diff: ndtiDiff,
    red_green_diff: rgDiff
  });
})).sort("system:time_start");

print("Monthly downstream-minus-upstream table", monthlyDiff);

var tssDiffChart = ui.Chart.feature.byFeature(monthlyDiff.filter(ee.Filter.neq("tss_diff", noDataValue)), "date", ["tss_diff"])
  .setChartType("ColumnChart")
  .setOptions({
    title: "Monthly TSS Difference (Downstream - Upstream)",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "Difference"},
    colors: ["#B71C1C"]
  });
print(tssDiffChart);

var ndtiDiffChart = ui.Chart.feature.byFeature(monthlyDiff.filter(ee.Filter.neq("ndti_diff", noDataValue)), "date", ["ndti_diff"])
  .setChartType("ColumnChart")
  .setOptions({
    title: "Monthly NDTI Difference (Downstream - Upstream)",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "Difference"},
    colors: ["#6A1B9A"]
  });
print(ndtiDiffChart);

var redGreenDiffChart = ui.Chart.feature.byFeature(monthlyDiff.filter(ee.Filter.neq("red_green_diff", noDataValue)), "date", ["red_green_diff"])
  .setChartType("ColumnChart")
  .setOptions({
    title: "Monthly Red/Green Difference (Downstream - Upstream)",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "Difference"},
    colors: ["#1B5E20"]
  });
print(redGreenDiffChart);
