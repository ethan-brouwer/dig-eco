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

print("Recent window", recentStart.format("YYYY-MM-dd"), recentEnd.format("YYYY-MM-dd"));
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
  var optical = img.select(["B2", "B3", "B4", "B8", "B11"]).multiply(0.0001);
  return img.addBands(optical, null, true);
}

function addMonth(img) {
  return img.set("month", ee.Date(img.get("system:time_start")).get("month"));
}

function addWaterBands(img) {
  var green = img.select("B3");
  var red = img.select("B4");
  var nir = img.select("B8");
  var swir1 = img.select("B11");

  var ndwi = img.normalizedDifference(["B3", "B8"]).rename("NDWI");
  var mndwi = img.normalizedDifference(["B3", "B11"]).rename("MNDWI");
  var ndti = img.normalizedDifference(["B4", "B3"]).rename("NDTI");
  var redGreen = red.divide(green).rename("RED_GREEN");
  var tssProxy = red.multiply(1000).rename("TSS_PROXY");

  return img.addBands([ndwi, mndwi, ndti, redGreen, tssProxy]);
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
print("Recent image count", recentCol.size());

var recent = recentCol.median().clip(aoi);
var recentObsCount = recentCol.select("B4").count().clip(aoi).rename("RECENT_OBS_COUNT");
var stableWater = waterMask(recent).selfMask().rename("STABLE_WATER");
var excludedFromWaterMask = permanentWater.and(stableWater.not()).selfMask().rename("EXCLUDED_FROM_WATER_MASK");

var waterDownstream = stableWater.and(downstreamCorridor).selfMask();
var waterUpstream = stableWater.and(upstreamCorridor).selfMask();

Map.addLayer(recent, {bands: ["B4", "B3", "B2"], min: 0.02, max: 0.30}, "Recent composite", false);
Map.addLayer(permanentWater.selfMask(), {palette: ["4FC3F7"]}, "Permanent water (JRC)", false);
Map.addLayer(stableWater, {palette: ["00BFFF"]}, "Current stable water", false);
Map.addLayer(excludedFromWaterMask, {palette: ["FF00FF"]}, "QA excluded from water mask", false);
Map.addLayer(waterDownstream, {palette: ["FF6D00"]}, "QA used downstream water pixels", false);
Map.addLayer(waterUpstream, {palette: ["2962FF"]}, "QA used upstream water pixels", false);
Map.addLayer(recentObsCount, {min: 0, max: 20, palette: ["2b2b2b", "f7f7f7", "00ff00"]}, "QA obs count", false);

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

  var stats = img.select(["TSS_PROXY", "NDTI"]).updateMask(maskImg).reduceRegion({
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
    ndti_mean: ee.Algorithms.If(enoughWater, stats.get("NDTI"), null)
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

var downstreamChartStats = downstreamStats.filter(ee.Filter.notNull(["tss_mean", "ndti_mean"]));

var tssDistanceChart = ui.Chart.feature.byFeature(downstreamChartStats, "label", ["tss_mean"])
  .setChartType("ColumnChart")
  .setOptions({
    title: "Past-Year TSS Proxy by Downstream Distance Bin",
    hAxis: {title: "Distance from outlet"},
    vAxis: {title: "TSS proxy (relative)"},
    colors: ["#D84315"]
  });
print(tssDistanceChart);

var ndtiDistanceChart = ui.Chart.feature.byFeature(downstreamChartStats, "label", ["ndti_mean"])
  .setChartType("ColumnChart")
  .setOptions({
    title: "Past-Year Water NDTI by Downstream Distance Bin",
    hAxis: {title: "Distance from outlet"},
    vAxis: {title: "NDTI"},
    colors: ["#8E24AA"]
  });
print(ndtiDistanceChart);

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

  var emptyFeature = ee.Feature(null, {
    "system:time_start": start.millis(),
    date: start.format("YYYY-MM"),
    img_count: count,
    up_tss: null,
    down_tss_0_500: null,
    down_tss_500_1500: null,
    down_tss_1500_3000: null
  });

  return ee.Feature(ee.Algorithms.If(count.gt(0), (function() {
    var img = ee.Image(col.median()).clip(aoi);
    var w = waterMask(img).selfMask();
    var upMask = w.and(upstreamCorridor).selfMask();

    function zoneTss(zoneFeature, corridorMask) {
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

      var mean = img.select("TSS_PROXY").updateMask(w).updateMask(corridorMask).reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: zoneGeom,
        scale: scaleMeters,
        bestEffort: true,
        maxPixels: 1e9
      }).get("TSS_PROXY");

      return ee.Algorithms.If(px.gte(minWaterPixelsPerZone), mean, null);
    }

    var downList = downstreamZones.toList(downstreamZones.size());
    return ee.Feature(null, {
      "system:time_start": start.millis(),
      date: start.format("YYYY-MM"),
      img_count: count,
      up_tss: zoneTss(upstreamZone, upMask),
      down_tss_0_500: zoneTss(ee.Feature(downList.get(0)), w.and(downstreamCorridor).selfMask()),
      down_tss_500_1500: zoneTss(ee.Feature(downList.get(1)), w.and(downstreamCorridor).selfMask()),
      down_tss_1500_3000: zoneTss(ee.Feature(downList.get(2)), w.and(downstreamCorridor).selfMask())
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

print("Simple downstream screening script ready. First verify the snapped outlet point, stream corridor, bin placement, and whether downstream monthly TSS stays above the upstream control.");
