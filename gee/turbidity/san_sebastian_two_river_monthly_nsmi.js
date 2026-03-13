/*
  FILE: gee/turbidity/san_sebastian_two_river_monthly_nsmi.js
  PURPOSE: Compare monthly NSMI mean for two river corridors using a water-only 10 m analysis mask.

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
var displayCorridorBufferMeters = 20;
var analysisCorridorBufferMeters = 10;
var analysisScaleMeters = 10;
var useWetSeasonFilter = false;
var wetSeasonMonths = [5, 6, 7, 8, 9, 10];
var ndwiThreshold = 0.0;
var mndwiThreshold = 0.0;
var nirWaterMax = 0.12;
var swirWaterMax = 0.10;
var minObsCount = 2;
var minWaterPixelsPerReach = 5;
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

var ssRiversDisplayCorridor = ssRiversGeom.buffer(displayCorridorBufferMeters);
var southRiverFullDisplayCorridor = southRiverFullGeom.buffer(displayCorridorBufferMeters);
var ssRiversAnalysisCorridor = ssRiversGeom.buffer(analysisCorridorBufferMeters);
var southRiverFullAnalysisCorridor = southRiverFullGeom.buffer(analysisCorridorBufferMeters);
var analysisCorridorsGeom = ssRiversAnalysisCorridor.union(southRiverFullAnalysisCorridor, 1);
var aoi = analysisCorridorsGeom
  .bounds()
  .buffer(500);

var comparisonReaches = ee.FeatureCollection([
  ee.Feature(ssRiversAnalysisCorridor, {
    reach_id: "SSrivers",
    reach_type: "reference"
  }),
  ee.Feature(southRiverFullAnalysisCorridor, {
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
Map.addLayer(ssRiversDisplayCorridor, {color: "26A69A"}, "SSrivers display corridor", false);
Map.addLayer(southRiverFullDisplayCorridor, {color: "E53935"}, "South_River_Full display corridor", false);
Map.addLayer(ssRiversAnalysisCorridor, {color: "00BCD4"}, "SSrivers analysis corridor", false);
Map.addLayer(southRiverFullAnalysisCorridor, {color: "FF7043"}, "South_River_Full analysis corridor", false);
Map.addLayer(aoi, {color: "FFFFFF"}, "AOI", false);

// ================= DATE WINDOW =================
var endDate = ee.Date(Date.now());
var startDate = ee.Date.fromYMD(analysisStartYear, 1, 1);

print("Analysis window", startDate.format("YYYY-MM-dd"), endDate.format("YYYY-MM-dd"));
print("Display corridor buffer (m)", displayCorridorBufferMeters);
print("Analysis corridor buffer (m)", analysisCorridorBufferMeters);
print("Analysis scale (m)", analysisScaleMeters);
print("Wet-season month filter", useWetSeasonFilter ? wetSeasonMonths : "OFF");
print("NDWI threshold", ndwiThreshold);
print("MNDWI threshold", mndwiThreshold);
print("NIR max (dark-water assist)", nirWaterMax);
print("SWIR1 max (dark-water assist)", swirWaterMax);
print("Min observations per month", minObsCount);
print("Min water pixels per reach", minWaterPixelsPerReach);

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

var analysisCorridorMask = ee.Image.constant(1)
  .clip(analysisCorridorsGeom)
  .selfMask()
  .rename("ANALYSIS_CORRIDOR_MASK");

function waterMask(img) {
  return img.select("NDWI").gt(ndwiThreshold)
    .or(img.select("MNDWI").gt(mndwiThreshold))
    .or(
      img.select("B8").lt(nirWaterMax)
        .and(img.select("B11").lt(swirWaterMax))
    )
    .updateMask(analysisCorridorMask)
    .selfMask()
    .rename("WATER_MASK");
}

// ================= MONTHLY NSMI COMPARISON =================
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

var monthlyNSMI = ee.FeatureCollection(monthlyStarts.map(function(mStart) {
  mStart = ee.Date(mStart);
  var mEnd = mStart.advance(1, "month");
  var col = s2.filterDate(mStart, mEnd);
  var count = col.size();

  var emptyFeatures = ee.FeatureCollection([
    ee.Feature(null, {
      "system:time_start": mStart.millis(),
      date: mStart.format("YYYY-MM"),
      year: ee.Number.parse(mStart.format("YYYY")),
      month_num: ee.Number.parse(mStart.format("M")),
      reach_id: "SSrivers",
      reach_type: "reference",
      image_count: count,
      water_px: 0,
      qa_flag: "no_images",
      nsmi_mean: noDataValue
    }),
    ee.Feature(null, {
      "system:time_start": mStart.millis(),
      date: mStart.format("YYYY-MM"),
      year: ee.Number.parse(mStart.format("YYYY")),
      month_num: ee.Number.parse(mStart.format("M")),
      reach_id: "South_River_Full",
      reach_type: "target",
      image_count: count,
      water_px: 0,
      qa_flag: "no_images",
      nsmi_mean: noDataValue
    })
  ]);

  return ee.FeatureCollection(ee.Algorithms.If(count.gt(0), (function() {
    var img = ee.Image(col.median()).clip(aoi);
    var monthObs = col.select("B4").count().clip(aoi);
    var monthWaterMask = waterMask(img)
      .updateMask(monthObs.gte(minObsCount))
      .selfMask();

    var ssWaterPx = safeNumber(monthWaterMask.reduceRegion({
      reducer: ee.Reducer.count(),
      geometry: ssRiversAnalysisCorridor,
      scale: analysisScaleMeters,
      bestEffort: true,
      maxPixels: 1e9
    }).values().get(0));

    var southWaterPx = safeNumber(monthWaterMask.reduceRegion({
      reducer: ee.Reducer.count(),
      geometry: southRiverFullAnalysisCorridor,
      scale: analysisScaleMeters,
      bestEffort: true,
      maxPixels: 1e9
    }).values().get(0));

    var waterNsmi = img.select("NSMI").updateMask(monthWaterMask);

    var ssStats = waterNsmi.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: ssRiversAnalysisCorridor,
      scale: analysisScaleMeters,
      bestEffort: true,
      maxPixels: 1e9
    });

    var southStats = waterNsmi.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: southRiverFullAnalysisCorridor,
      scale: analysisScaleMeters,
      bestEffort: true,
      maxPixels: 1e9
    });

    return ee.FeatureCollection([
      ee.Feature(null, {
        "system:time_start": mStart.millis(),
        date: mStart.format("YYYY-MM"),
        year: ee.Number.parse(mStart.format("YYYY")),
        month_num: ee.Number.parse(mStart.format("M")),
        reach_id: "SSrivers",
        reach_type: "reference",
        image_count: count,
        water_px: ssWaterPx,
        qa_flag: ee.Algorithms.If(ssWaterPx.gte(minWaterPixelsPerReach), "water_mask", "low_water_px"),
        nsmi_mean: ee.Algorithms.If(
          ssWaterPx.gte(minWaterPixelsPerReach),
          safeNumber(ssStats.get("NSMI")),
          noDataValue
        )
      }),
      ee.Feature(null, {
        "system:time_start": mStart.millis(),
        date: mStart.format("YYYY-MM"),
        year: ee.Number.parse(mStart.format("YYYY")),
        month_num: ee.Number.parse(mStart.format("M")),
        reach_id: "South_River_Full",
        reach_type: "target",
        image_count: count,
        water_px: southWaterPx,
        qa_flag: ee.Algorithms.If(southWaterPx.gte(minWaterPixelsPerReach), "water_mask", "low_water_px"),
        nsmi_mean: ee.Algorithms.If(
          southWaterPx.gte(minWaterPixelsPerReach),
          safeNumber(southStats.get("NSMI")),
          noDataValue
        )
      })
    ]);
  })(), emptyFeatures));
})).flatten().sort("system:time_start");

print("Monthly NSMI table", monthlyNSMI);
print("Comparison reaches", comparisonReaches);

// ================= CHARTS =================
var monthlyNSMIChartFc = monthlyNSMI.map(function(f) {
  var meanVal = f.get("nsmi_mean");
  var waterPx = f.get("water_px");
  return f.set({
    nsmi_mean_chart: ee.Algorithms.If(ee.Algorithms.IsEqual(meanVal, noDataValue), null, meanVal),
    water_px_chart: ee.Algorithms.If(ee.Number(waterPx).gt(0), waterPx, null)
  });
});

var nsmiMeanChart = ui.Chart.feature.groups(
  monthlyNSMIChartFc.filter(ee.Filter.notNull(["nsmi_mean_chart"])),
  "date",
  "nsmi_mean_chart",
  "reach_id"
).setChartType("LineChart")
  .setOptions({
    title: "Monthly NSMI Mean",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "NSMI mean"},
    lineWidth: 2,
    pointSize: 3,
    colors: ["#26A69A", "#E53935"]
  });
print(nsmiMeanChart);

var waterPxChart = ui.Chart.feature.groups(
  monthlyNSMIChartFc.filter(ee.Filter.notNull(["water_px_chart"])),
  "date",
  "water_px_chart",
  "reach_id"
).setChartType("LineChart")
  .setOptions({
    title: "Monthly Water Pixels Used",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "Water pixels at 10 m"},
    lineWidth: 2,
    pointSize: 3,
    colors: ["#26A69A", "#E53935"]
  });
print(waterPxChart);
