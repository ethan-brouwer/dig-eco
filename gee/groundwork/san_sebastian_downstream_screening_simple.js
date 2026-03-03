/*
  FILE: gee/groundwork/san_sebastian_downstream_screening_simple.js
  PURPOSE: Minimal past-year downstream screening workflow for San Sebastian Mine
           using HydroSHEDS stream extraction plus Sentinel-2 water proxies.

  DESIGN GOAL
  - Keep the workflow simple enough to sanity-check in the GEE Code Editor.
  - Use only a few screening layers and charts before building a heavier workflow.
  - Detect possible downstream sediment/discoloration signals, not mercury directly.

  IMPORTANT LIMITATIONS
  - This script does NOT measure dissolved mercury.
  - The stream path is HydroSHEDS-based, but the downstream side is still constrained
    by a user-set azimuth for a lightweight first-pass implementation.
  - Outputs are screening proxies only and should be validated with field data.
*/

// ================= USER SETTINGS =================
var siteName = "San Sebastian Mine (MRDS)";
var siteLon = -87.93002;
var siteLat = 13.6509;

// Keep AOI tight for fast testing.
var aoiRadiusMeters = 5000;
var recentMonths = 12;
var baselineYears = 3;
var cloudMax = 60;
var seasonMonths = [11, 12, 1, 2, 3, 4];
var useSeasonFilter = false;

// Hydro / stream settings.
var downstreamAzimuthDeg = 110;    // first-pass guide; adjust after visual inspection
var streamAccThreshold = 250;      // lower = more channels, higher = main stem only
var streamBufferMeters = 60;
var outletSnapMeters = 300;
var distanceBinEdgesMeters = [0, 500, 1500, 3000];
var upstreamControlRangeMeters = [-1500, -500];

// QA / scale settings.
var scaleMeters = 30;              // use 30 m for stable HydroSHEDS + S2 summary
var minWaterPixelsPerZone = 5;

// ================= MAP SETUP =================
Map.setOptions("SATELLITE");

var sitePoint = ee.Geometry.Point([siteLon, siteLat]);
var aoi = sitePoint.buffer(aoiRadiusMeters);
Map.centerObject(sitePoint, 13);

Map.addLayer(sitePoint, {color: "FF0000"}, siteName, true);
Map.addLayer(aoi, {color: "FFFFFF"}, "AOI (5 km)", false);

// ================= DATE WINDOW =================
var now = ee.Date(Date.now());
var recentStart = now.advance(-recentMonths, "month");
var recentEnd = now;
var baselineStart = recentStart.advance(-baselineYears, "year");
var baselineEnd = recentStart;

print("Recent window", recentStart.format("YYYY-MM-dd"), recentEnd.format("YYYY-MM-dd"));
print("Baseline window", baselineStart.format("YYYY-MM-dd"), baselineEnd.format("YYYY-MM-dd"));
print("Distance bins (m)", distanceBinEdgesMeters);

// ================= SENTINEL-2 HELPERS =================
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
  var optical = img.select(["B2", "B3", "B4", "B5", "B8", "B11"]).multiply(0.0001);
  return img.addBands(optical, null, true);
}

function addMonth(img) {
  return img.set("month", ee.Date(img.get("system:time_start")).get("month"));
}

function addWaterBands(img) {
  var green = img.select("B3");
  var red = img.select("B4");
  var redEdge1 = img.select("B5");
  var nir = img.select("B8");
  var swir1 = img.select("B11");

  var ndwi = img.normalizedDifference(["B3", "B8"]).rename("NDWI");
  var mndwi = img.normalizedDifference(["B3", "B11"]).rename("MNDWI");
  var ndti = img.normalizedDifference(["B4", "B3"]).rename("NDTI");
  var redGreen = red.divide(green).rename("RED_GREEN");
  var tssProxy = red.multiply(1000).rename("TSS_PROXY");
  var redEdgeTurb = redEdge1.rename("RED_EDGE_TURB");

  return img.addBands([ndwi, mndwi, ndti, redGreen, tssProxy, redEdgeTurb]);
}

function buildS2Collection(start, end) {
  var col = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
    .filterBounds(aoi)
    .filterDate(start, end)
    .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", cloudMax))
    .map(maskS2)
    .map(scaleS2)
    .map(addMonth);

  if (useSeasonFilter) {
    col = col.filter(ee.Filter.inList("month", seasonMonths));
  }

  return col.map(addWaterBands);
}

function waterMask(img) {
  var opticalWater = img.select("NDWI").gt(0.12).or(img.select("MNDWI").gt(0.15));
  return opticalWater.and(permanentWater);
}

// ================= HYDROSHEDS STREAM SCREEN =================
var hydroAcc = ee.Image("WWF/HydroSHEDS/15ACC").clip(aoi);
var streamMask = hydroAcc.gte(streamAccThreshold).rename("STREAM_MASK");

var azimuthRad = downstreamAzimuthDeg * Math.PI / 180;
var ux = Math.sin(azimuthRad);
var uy = Math.cos(azimuthRad);

var lonLat = ee.Image.pixelLonLat();
var signedAlong = lonLat.select("longitude").subtract(siteLon).multiply(ux)
  .add(lonLat.select("latitude").subtract(siteLat).multiply(uy))
  .rename("SIGNED_ALONG");

var snapCircle = sitePoint.buffer(outletSnapMeters);
var snappedPointFc = streamMask.selfMask().reduceToVectors({
  geometry: snapCircle,
  scale: scaleMeters,
  geometryType: "centroid",
  maxPixels: 1e9
}).map(function(f) {
  return f.set("dist_m", f.geometry().distance(sitePoint, 1));
});

var snappedPointCount = snappedPointFc.size();
var snappedPoint = ee.Feature(snappedPointFc.sort("dist_m").first());
var snappedGeom = ee.Geometry(ee.Algorithms.If(
  snappedPointCount.gt(0),
  snappedPoint.geometry(),
  sitePoint
));

var sourceImage = ee.Image().toByte().paint(snappedGeom, 1).selfMask();
var streamCost = ee.Image.constant(1).updateMask(streamMask);
var streamDistance = streamCost.cumulativeCost(sourceImage, aoiRadiusMeters)
  .rename("STREAM_DISTANCE_M");

var downstreamGuide = signedAlong.gte(0);
var upstreamGuide = signedAlong.lt(0);

var downstreamStream = streamMask.and(downstreamGuide).and(streamDistance.lte(aoiRadiusMeters));
var upstreamStream = streamMask.and(upstreamGuide).and(streamDistance.lte(Math.abs(upstreamControlRangeMeters[0])));

var downstreamCorridor = downstreamStream.focal_max(streamBufferMeters, "circle", "meters")
  .selfMask()
  .rename("DOWNSTREAM_CORRIDOR");

var upstreamCorridor = upstreamStream.focal_max(streamBufferMeters, "circle", "meters")
  .selfMask()
  .rename("UPSTREAM_CORRIDOR");

Map.addLayer(streamMask.selfMask(), {palette: ["00FFFF"]}, "HydroSHEDS stream mask", false);
Map.addLayer(snappedGeom, {color: "FFFF00"}, "Snapped outlet point", true);
Map.addLayer(downstreamCorridor, {palette: ["FF6D00"]}, "Downstream corridor", true);
Map.addLayer(upstreamCorridor, {palette: ["2962FF"]}, "Upstream control corridor", true);

// ================= WATER COMPOSITE =================
var jrcWater = ee.Image("JRC/GSW1_4/GlobalSurfaceWater").select("occurrence").clip(aoi);
var permanentWater = jrcWater.gte(50).rename("PERMANENT_WATER");

var recentCol = buildS2Collection(recentStart, recentEnd);
var baselineCol = buildS2Collection(baselineStart, baselineEnd);
print("Recent image count", recentCol.size());
print("Baseline image count", baselineCol.size());

var recent = recentCol.median().clip(aoi);
var baseline = baselineCol.median().clip(aoi);
var recentObsCount = recentCol.select("B4").count().clip(aoi).rename("RECENT_OBS_COUNT");
var baselineObsCount = baselineCol.select("B4").count().clip(aoi).rename("BASELINE_OBS_COUNT");
var recentWater = waterMask(recent).selfMask();
var baselineWater = waterMask(baseline).selfMask();
var stableWater = recentWater.and(baselineWater).selfMask().rename("STABLE_WATER");
var excludedFromWaterMask = permanentWater.and(stableWater.not()).selfMask().rename("EXCLUDED_FROM_WATER_MASK");

var tssAnomaly = recent.select("TSS_PROXY").subtract(baseline.select("TSS_PROXY"))
  .updateMask(stableWater)
  .rename("TSS_ANOMALY");
var ndtiAnomaly = recent.select("NDTI").subtract(baseline.select("NDTI"))
  .updateMask(stableWater)
  .rename("NDTI_ANOMALY");
var redEdgeAnomaly = recent.select("RED_EDGE_TURB").subtract(baseline.select("RED_EDGE_TURB"))
  .updateMask(stableWater)
  .rename("RED_EDGE_TURB_ANOMALY");

var waterDownstream = stableWater.and(downstreamCorridor).selfMask();
var waterUpstream = stableWater.and(upstreamCorridor).selfMask();

Map.addLayer(recent, {bands: ["B4", "B3", "B2"], min: 0.02, max: 0.30}, "Recent composite", false);
Map.addLayer(permanentWater.selfMask(), {palette: ["4FC3F7"]}, "Permanent water (JRC)", false);
Map.addLayer(stableWater, {palette: ["00BFFF"]}, "Current stable water", false);
Map.addLayer(excludedFromWaterMask, {palette: ["FF00FF"]}, "QA excluded from water mask", false);
Map.addLayer(waterDownstream, {palette: ["FF6D00"]}, "QA used downstream water pixels", false);
Map.addLayer(waterUpstream, {palette: ["2962FF"]}, "QA used upstream water pixels", false);
Map.addLayer(recentObsCount, {min: 0, max: 20, palette: ["2b2b2b", "f7f7f7", "00ff00"]}, "QA obs count", false);
Map.addLayer(baselineObsCount, {min: 0, max: 40, palette: ["2b2b2b", "f7f7f7", "00ff00"]}, "QA baseline obs count", false);
Map.addLayer(tssAnomaly, {min: -40, max: 40, palette: ["2166AC", "F7F7F7", "B2182B"]}, "TSS anomaly vs baseline", false);
Map.addLayer(ndtiAnomaly, {min: -0.15, max: 0.15, palette: ["2166AC", "F7F7F7", "B2182B"]}, "NDTI anomaly vs baseline", false);
Map.addLayer(redEdgeAnomaly, {min: -0.03, max: 0.03, palette: ["2166AC", "F7F7F7", "B2182B"]}, "Red-edge turbidity anomaly vs baseline", false);

// ================= DISTANCE BINS =================
var distanceBreaks = ee.List(distanceBinEdgesMeters);
var downstreamZones = ee.FeatureCollection(ee.List.sequence(0, distanceBreaks.length().subtract(2)).map(function(i) {
  i = ee.Number(i);
  var inner = ee.Number(distanceBreaks.get(i));
  var outer = ee.Number(distanceBreaks.get(i.add(1)));
  var zoneMask = streamDistance.gte(inner).and(streamDistance.lt(outer)).and(downstreamCorridor);
  var zoneGeom = zoneMask.selfMask().reduceToVectors({
    geometry: aoi,
    scale: scaleMeters,
    geometryType: "polygon",
    maxPixels: 1e9
  }).geometry();

  return ee.Feature(zoneGeom, {
    zone_id: ee.String("down_").cat(inner.format()).cat("_").cat(outer.format()),
    zone_type: "downstream",
    inner_m: inner,
    outer_m: outer,
    label: inner.format().cat("-").cat(outer.format()).cat(" m")
  });
}));

var upstreamZoneMask = streamDistance.gte(Math.abs(upstreamControlRangeMeters[1]))
  .and(streamDistance.lt(Math.abs(upstreamControlRangeMeters[0])))
  .and(upstreamCorridor);

var upstreamZone = ee.Feature(
  upstreamZoneMask.selfMask().reduceToVectors({
    geometry: aoi,
    scale: scaleMeters,
    geometryType: "polygon",
    maxPixels: 1e9
  }).geometry(),
  {
    zone_id: "upstream_control",
    zone_type: "upstream",
    inner_m: upstreamControlRangeMeters[0],
    outer_m: upstreamControlRangeMeters[1],
    label: "Upstream control"
  }
);

Map.addLayer(downstreamZones.style({
  color: "FFA500",
  fillColor: "00000000",
  width: 2
}), {}, "Downstream bins", true);

Map.addLayer(ee.FeatureCollection([upstreamZone]).style({
  color: "2962FF",
  fillColor: "00000000",
  width: 2
}), {}, "Upstream control zone", true);

// Bin-specific QA layers so the user can verify exactly which pixels feed each summary.
var downstreamZoneList = downstreamZones.toList(downstreamZones.size());
var downZone0 = ee.Feature(downstreamZoneList.get(0));
var downZone1 = ee.Feature(downstreamZoneList.get(1));
var downZone2 = ee.Feature(downstreamZoneList.get(2));

var downZone0Pixels = waterDownstream.clip(downZone0.geometry()).selfMask().rename("DOWN_0_500_PX");
var downZone1Pixels = waterDownstream.clip(downZone1.geometry()).selfMask().rename("DOWN_500_1500_PX");
var downZone2Pixels = waterDownstream.clip(downZone2.geometry()).selfMask().rename("DOWN_1500_3000_PX");
var upstreamZonePixels = waterUpstream.clip(upstreamZone.geometry()).selfMask().rename("UPSTREAM_CONTROL_PX");

Map.addLayer(downZone0Pixels, {palette: ["FF0000"]}, "QA bin pixels 0-500 m", false);
Map.addLayer(downZone1Pixels, {palette: ["FF9800"]}, "QA bin pixels 500-1500 m", false);
Map.addLayer(downZone2Pixels, {palette: ["FFD54F"]}, "QA bin pixels 1500-3000 m", false);
Map.addLayer(upstreamZonePixels, {palette: ["2962FF"]}, "QA upstream control pixels", false);

// ================= ZONAL SUMMARIES =================
function summarizeZone(feature, img, maskImg) {
  var geom = feature.geometry();
  var zoneMask = ee.Image.constant(1).updateMask(maskImg).clip(geom);
  var waterDict = zoneMask.reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: geom,
    scale: scaleMeters,
    bestEffort: true,
    maxPixels: 1e9
  });

  var waterPxRaw = waterDict.values().get(0);
  var waterPx = ee.Number(ee.Algorithms.If(
    ee.Algorithms.IsEqual(waterPxRaw, null),
    0,
    waterPxRaw
  ));

  var stats = img.select(["TSS_PROXY", "NDTI", "RED_EDGE_TURB"]).updateMask(maskImg).reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: geom,
    scale: scaleMeters,
    bestEffort: true,
    maxPixels: 1e9
  });

  var enoughWater = waterPx.gte(minWaterPixelsPerZone);

  return feature.set({
    water_px: waterPx,
    tss_mean: ee.Algorithms.If(enoughWater, stats.get("TSS_PROXY"), null),
    ndti_mean: ee.Algorithms.If(enoughWater, stats.get("NDTI"), null),
    red_edge_mean: ee.Algorithms.If(enoughWater, stats.get("RED_EDGE_TURB"), null)
  });
}

var downstreamStats = downstreamZones.map(function(f) {
  return summarizeZone(f, recent, waterDownstream);
});

var upstreamStats = ee.FeatureCollection([
  summarizeZone(upstreamZone, recent, waterUpstream)
]);

var profileStats = upstreamStats.merge(downstreamStats);
print("Downstream profile stats", profileStats);

var downstreamTssChartStats = downstreamStats.filter(ee.Filter.notNull(["tss_mean"]));
var downstreamNdtiChartStats = downstreamStats.filter(ee.Filter.notNull(["ndti_mean"]));
var downstreamRedEdgeChartStats = downstreamStats.filter(ee.Filter.notNull(["red_edge_mean"]));

var tssDistanceChart = ui.Chart.feature.byFeature(downstreamTssChartStats, "label", ["tss_mean"])
  .setChartType("ColumnChart")
  .setOptions({
    title: "Past-Year TSS Proxy by Downstream Distance Bin",
    hAxis: {title: "Distance from outlet"},
    vAxis: {title: "TSS proxy (relative)"},
    colors: ["#D84315"]
  });
print(tssDistanceChart);

var ndtiDistanceChart = ui.Chart.feature.byFeature(downstreamNdtiChartStats, "label", ["ndti_mean"])
  .setChartType("ColumnChart")
  .setOptions({
    title: "Past-Year Water NDTI by Downstream Distance Bin",
    hAxis: {title: "Distance from outlet"},
    vAxis: {title: "NDTI"},
    colors: ["#8E24AA"]
  });
print(ndtiDistanceChart);

var redEdgeDistanceChart = ui.Chart.feature.byFeature(downstreamRedEdgeChartStats, "label", ["red_edge_mean"])
  .setChartType("ColumnChart")
  .setOptions({
    title: "Past-Year Red-Edge Turbidity Proxy by Downstream Distance Bin",
    hAxis: {title: "Distance from outlet"},
    vAxis: {title: "Red-edge reflectance (relative)"},
    colors: ["#2E7D32"]
  });
print(redEdgeDistanceChart);

// ================= SIMPLE MONTHLY CONTROL CHART =================
function monthStartList(start, end) {
  var monthCount = ee.Number(end.difference(start, "month")).floor();
  return ee.List.sequence(0, monthCount.subtract(1)).map(function(m) {
    return start.advance(ee.Number(m), "month");
  });
}

var monthlyStarts = monthStartList(ee.Date(recentStart.format("YYYY-MM-01")), ee.Date(recentEnd.format("YYYY-MM-01")).advance(1, "month"));

var monthlyStats = ee.FeatureCollection(monthlyStarts.map(function(start) {
  start = ee.Date(start);
  var end = start.advance(1, "month");
  var col = buildS2Collection(start, end);
  var count = col.size();
  var rainMm = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY")
    .filterBounds(aoi)
    .filterDate(start, end)
    .sum()
    .reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: aoi,
      scale: 5000,
      bestEffort: true,
      maxPixels: 1e9
    }).get("precipitation");

  var emptyFeature = ee.Feature(null, {
    "system:time_start": start.millis(),
    date: start.format("YYYY-MM"),
    img_count: count,
    rain_mm: rainMm,
    up_tss: null,
    down_tss_0_500: null,
    down_tss_500_1500: null,
    down_tss_1500_3000: null,
    up_red_edge: null,
    down_red_edge_0_500: null,
    down_minus_up_tss_0_500: null
  });

  return ee.Feature(ee.Algorithms.If(count.gt(0), (function() {
    var img = ee.Image(col.median()).clip(aoi);
    var w = waterMask(img).selfMask();
    var upMask = w.and(upstreamCorridor).selfMask();

    function zoneMean(zoneFeature, corridorMask, bandName) {
      var zoneGeom = zoneFeature.geometry();
      var zoneMask = ee.Image.constant(1).updateMask(w).updateMask(corridorMask).clip(zoneGeom);
      var pxDict = zoneMask.reduceRegion({
        reducer: ee.Reducer.count(),
        geometry: zoneGeom,
        scale: scaleMeters,
        bestEffort: true,
        maxPixels: 1e9
      });
      var pxRaw = pxDict.values().get(0);
      var px = ee.Number(ee.Algorithms.If(
        ee.Algorithms.IsEqual(pxRaw, null),
        0,
        pxRaw
      ));

      var mean = img.select(bandName).updateMask(w).updateMask(corridorMask).reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: zoneGeom,
        scale: scaleMeters,
        bestEffort: true,
        maxPixels: 1e9
      }).get(bandName);

      return ee.Algorithms.If(px.gte(minWaterPixelsPerZone), mean, null);
    }

    var downList = downstreamZones.toList(downstreamZones.size());
    var upTss = zoneMean(upstreamZone, upMask, "TSS_PROXY");
    var downTss0 = zoneMean(ee.Feature(downList.get(0)), w.and(downstreamCorridor).selfMask(), "TSS_PROXY");
    var upRedEdge = zoneMean(upstreamZone, upMask, "RED_EDGE_TURB");
    var downRedEdge0 = zoneMean(ee.Feature(downList.get(0)), w.and(downstreamCorridor).selfMask(), "RED_EDGE_TURB");

    return ee.Feature(null, {
      "system:time_start": start.millis(),
      date: start.format("YYYY-MM"),
      img_count: count,
      rain_mm: rainMm,
      up_tss: upTss,
      down_tss_0_500: downTss0,
      down_tss_500_1500: zoneMean(ee.Feature(downList.get(1)), w.and(downstreamCorridor).selfMask(), "TSS_PROXY"),
      down_tss_1500_3000: zoneMean(ee.Feature(downList.get(2)), w.and(downstreamCorridor).selfMask(), "TSS_PROXY"),
      up_red_edge: upRedEdge,
      down_red_edge_0_500: downRedEdge0,
      down_minus_up_tss_0_500: ee.Algorithms.If(
        ee.Algorithms.IsEqual(upTss, null),
        null,
        ee.Algorithms.If(
          ee.Algorithms.IsEqual(downTss0, null),
          null,
          ee.Number(downTss0).subtract(ee.Number(upTss))
        )
      )
    });
  })(), emptyFeature));
})).sort("system:time_start");

print("Monthly stats table", monthlyStats);

var monthlyTssChart = ui.Chart.feature.byFeature(
  monthlyStats,
  "date",
  ["up_tss", "down_tss_0_500", "down_tss_500_1500", "down_tss_1500_3000"]
).setChartType("LineChart")
  .setOptions({
    title: "Monthly TSS Proxy: Upstream Control vs Downstream Bins",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "TSS proxy (relative)"},
    lineWidth: 2,
    pointSize: 3,
    series: {
      0: {color: "#1565C0"},
      1: {color: "#E53935"},
      2: {color: "#FB8C00"},
      3: {color: "#6D4C41"}
    }
  });
print(monthlyTssChart);

var monthlyTssDiffChart = ui.Chart.feature.byFeature(monthlyStats, "date", ["down_minus_up_tss_0_500"])
  .setChartType("ColumnChart")
  .setOptions({
    title: "Monthly TSS Proxy Difference (0-500 m downstream minus upstream)",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "Difference (relative)"},
    colors: ["#B71C1C"]
  });
print(monthlyTssDiffChart);

var monthlyRedEdgeChart = ui.Chart.feature.byFeature(
  monthlyStats,
  "date",
  ["up_red_edge", "down_red_edge_0_500"]
).setChartType("LineChart")
  .setOptions({
    title: "Monthly Red-Edge Turbidity Proxy: Upstream vs 0-500 m Downstream",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "Red-edge reflectance (relative)"},
    lineWidth: 2,
    pointSize: 3,
    series: {
      0: {color: "#1565C0"},
      1: {color: "#2E7D32"}
    }
  });
print(monthlyRedEdgeChart);

var monthlyRainChart = ui.Chart.feature.byFeature(monthlyStats, "date", ["rain_mm"])
  .setChartType("ColumnChart")
  .setOptions({
    title: "Monthly CHIRPS Rainfall (AOI mean, mm)",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "Rainfall (mm)"},
    colors: ["#26A69A"]
  });
print(monthlyRainChart);

var monthlyQaChart = ui.Chart.feature.byFeature(monthlyStats, "date", ["img_count"])
  .setChartType("ColumnChart")
  .setOptions({
    title: "QA: monthly Sentinel-2 image count",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "Image count"},
    colors: ["#546E7A"]
  });
print(monthlyQaChart);

// ================= OPTIONAL EXPORTS (UNCOMMENT IF NEEDED) =================
// Export.table.toDrive({
//   collection: profileStats,
//   description: "san_sebastian_downstream_profile_stats_simple",
//   folder: "GEE_exports",
//   fileNamePrefix: "san_sebastian_downstream_profile_stats_simple",
//   fileFormat: "CSV"
// });
//
// Export.table.toDrive({
//   collection: monthlyStats,
//   description: "san_sebastian_downstream_monthly_tss_simple",
//   folder: "GEE_exports",
//   fileNamePrefix: "san_sebastian_downstream_monthly_tss_simple",
//   fileFormat: "CSV"
// });

print("Simple downstream screening script ready. Verify the snapped outlet point, stream corridor, bin placement, QA pixel layers, baseline anomaly layers, and whether downstream TSS/red-edge signals rise above the upstream control during wetter months.");
