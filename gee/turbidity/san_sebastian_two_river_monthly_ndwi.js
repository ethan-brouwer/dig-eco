/*
  FILE: gee/turbidity/san_sebastian_two_river_monthly_ndwi.js
  PURPOSE: Compare monthly NDWI mean and median for two buffered river corridors.

  GEE IMPORTS REQUIRED
  - impact_point
  - SSrivers
  - South_River_Full
*/

// ================= USER SETTINGS =================
var siteName = "San Sebastian Mine (MRDS)";
var siteLon = -87.928972;
var siteLat = 13.646692;

var analysisStartYear = 2017;
var cloudMax = 60;
var corridorBufferMeters = 20;
var useWetSeasonFilter = false;
var wetSeasonMonths = [5, 6, 7, 8, 9, 10];
var noDataValue = -9999;

// ================= GEOMETRY HELPERS =================
var sitePoint = ee.Geometry.Point([siteLon, siteLat]);

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
var ssRiversGeom = toGeometry(SSrivers);
var southRiverFullGeom = toGeometry(South_River_Full);

var ssRiversCorridor = ssRiversGeom.buffer(corridorBufferMeters);
var southRiverFullCorridor = southRiverFullGeom.buffer(corridorBufferMeters);

// Union is only used to limit processing extent.
var riverMaskGeom = ssRiversCorridor.union(southRiverFullCorridor, 1);
var aoi = riverMaskGeom.bounds().buffer(500);

var comparisonReaches = ee.FeatureCollection([
  ee.Feature(ssRiversCorridor, {
    reach_id: "SSrivers",
    reach_type: "reference"
  }),
  ee.Feature(southRiverFullCorridor, {
    reach_id: "South_River_Full",
    reach_type: "target"
  })
]);

// ================= MAP SETUP =================
Map.setOptions("SATELLITE");
Map.centerObject(impactGeom, 12);
Map.addLayer(sitePoint, {color: "FF0000"}, siteName, false);
Map.addLayer(impactGeom, {color: "FFFF00"}, "impact_point", true);
Map.addLayer(ssRiversGeom, {color: "26A69A"}, "SSrivers line", true);
Map.addLayer(southRiverFullGeom, {color: "E53935"}, "South_River_Full line", true);
Map.addLayer(ssRiversCorridor, {color: "26A69A"}, "SSrivers corridor", false);
Map.addLayer(southRiverFullCorridor, {color: "E53935"}, "South_River_Full corridor", false);
Map.addLayer(aoi, {color: "FFFFFF"}, "AOI", false);

// ================= DATE WINDOW =================
var endDate = ee.Date(Date.now());
var startDate = ee.Date.fromYMD(analysisStartYear, 1, 1);

print("Analysis window", startDate.format("YYYY-MM-dd"), endDate.format("YYYY-MM-dd"));
print("Corridor buffer (m)", corridorBufferMeters);
print("Wet-season month filter", useWetSeasonFilter ? wetSeasonMonths : "OFF");

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
  img = img.clip(riverMaskGeom);
  var ndwi = img.normalizedDifference(["B3", "B8"]).rename("NDWI");
  return img.addBands([ndwi]);
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

// ================= MONTHLY NDWI COMPARISON =================
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

var monthlyNDWI = ee.FeatureCollection(monthlyStarts.map(function(mStart) {
  mStart = ee.Date(mStart);
  var mEnd = mStart.advance(1, "month");
  var col = s2.filterDate(mStart, mEnd);
  var count = col.size();

  return comparisonReaches.map(function(reach) {
    reach = ee.Feature(reach);
    var rGeom = reach.geometry();

    var baseFeature = ee.Feature(null, {
      "system:time_start": mStart.millis(),
      date: mStart.format("YYYY-MM"),
      year: ee.Number.parse(mStart.format("YYYY")),
      month_num: ee.Number.parse(mStart.format("M")),
      reach_id: reach.get("reach_id"),
      reach_type: reach.get("reach_type"),
      image_count: count,
      ndwi_mean: noDataValue,
      ndwi_median: noDataValue
    });

    return ee.Feature(ee.Algorithms.If(count.gt(0), (function() {
      var img = ee.Image(col.median()).clip(rGeom);
      var stats = img.select("NDWI").reduceRegion({
        reducer: ee.Reducer.mean().combine({
          reducer2: ee.Reducer.median(),
          sharedInputs: true
        }),
        geometry: rGeom,
        scale: 10,
        bestEffort: true,
        maxPixels: 1e9
      });

      return baseFeature.set({
        ndwi_mean: safeNumber(stats.get("NDWI_mean")),
        ndwi_median: safeNumber(stats.get("NDWI_median"))
      });
    })(), baseFeature));
  });
})).flatten().sort("system:time_start");

print("Monthly NDWI table", monthlyNDWI);

// ================= CHARTS =================
var monthlyNDWIChartFc = monthlyNDWI.map(function(f) {
  var meanVal = f.get("ndwi_mean");
  var medianVal = f.get("ndwi_median");
  return f.set({
    ndwi_mean_chart: ee.Algorithms.If(ee.Algorithms.IsEqual(meanVal, noDataValue), null, meanVal),
    ndwi_median_chart: ee.Algorithms.If(ee.Algorithms.IsEqual(medianVal, noDataValue), null, medianVal)
  });
});

var ndwiMeanChart = ui.Chart.feature.groups(
  monthlyNDWIChartFc.filter(ee.Filter.notNull(["ndwi_mean_chart"])),
  "date",
  "ndwi_mean_chart",
  "reach_id"
).setChartType("LineChart")
  .setOptions({
    title: "Monthly NDWI Mean",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "NDWI mean"},
    lineWidth: 2,
    pointSize: 3,
    colors: ["#26A69A", "#E53935"]
  });
print(ndwiMeanChart);

var ndwiMedianChart = ui.Chart.feature.groups(
  monthlyNDWIChartFc.filter(ee.Filter.notNull(["ndwi_median_chart"])),
  "date",
  "ndwi_median_chart",
  "reach_id"
).setChartType("LineChart")
  .setOptions({
    title: "Monthly NDWI Median",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "NDWI median"},
    lineWidth: 2,
    pointSize: 3,
    colors: ["#26A69A", "#E53935"]
  });
print(ndwiMedianChart);
