/*
  FILE: gee/groundwork/san_sebastian_ndvi_progressive_buffers.js
  PURPOSE: Lightweight NDVI site characterization around a selected site using
           progressively larger buffers/rings in Google Earth Engine.

  DESIGN GOAL
  - Fast(er) and modular: NDVI only (no anomaly scoring, no multi-index stacks).
  - Reusable for other mine sites by changing site coordinates and buffer distances.
  - Produces map layers + buffer/ring NDVI summary tables + simple charts.
*/

// ================= USER SETTINGS =================
var siteName = "San Sebastian Mine (MRDS)";
var siteLon = -87.93002;
var siteLat = 13.6509;

// Progressive distances from site center (meters).
// Script computes BOTH cumulative buffers (0-distance) and ring annuli.
var bufferDistancesMeters = [250, 500, 1000, 2000, 3000];

// Temporal / QA settings
// Use explicit dates for reproducibility across runs/sites.
// (Change these values directly when reusing the script.)
var startDate = "2025-01-01";
var endDate = "2026-01-01";   // exclusive end date in GEE filterDate

var cloudMax = 60;
var seasonMonths = [11, 12, 1, 2, 3, 4]; // drier season often cleaner observations
var useSeasonFilter = true;

// Analysis settings
var scaleMeters = 10;   // Sentinel-2 native for NDVI
var minObsCount = 3;    // stronger default QA: require at least 3 observations
var lowNdviThreshold = 0.2; // for "low vegetation" fraction by zone
var excludeBuiltUpAndWater = true; // recommended for mine-impact NDVI interpretation

// Metadata fields (kept simple so all future index scripts can match this structure)
var dataSource = "COPERNICUS/S2_SR_HARMONIZED";
var metricName = "NDVI";
var compositeMethod = "median";

// ================= MAP SETUP =================
Map.setOptions("SATELLITE");

var sitePoint = ee.Geometry.Point([siteLon, siteLat]);
var maxRadius = ee.Number(ee.List(bufferDistancesMeters).reduce(ee.Reducer.max()));
var aoi = sitePoint.buffer(maxRadius);

Map.centerObject(sitePoint, 13);
Map.addLayer(sitePoint, {color: "FF0000"}, siteName, true);
Map.addLayer(aoi, {color: "FFFFFF"}, "AOI (max buffer)", true);

// ================= DATE WINDOW =================
var recentStart = ee.Date(startDate);
var recentEnd = ee.Date(endDate);

print("Analysis window", recentStart.format("YYYY-MM-dd"), recentEnd.format("YYYY-MM-dd"));
print("Buffer distances (m)", bufferDistancesMeters);
print("Data source", dataSource);
print("Metric", metricName);
print("Composite method", compositeMethod);
print("Season filter (months)", useSeasonFilter ? seasonMonths : "OFF");

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
  var optical = img.select(["B4", "B8"]).multiply(0.0001);
  return img.addBands(optical, null, true);
}

function addMonth(img) {
  return img.set("month", ee.Date(img.get("system:time_start")).get("month"));
}

function addNdvi(img) {
  return img.addBands(
    img.normalizedDifference(["B8", "B4"]).rename("NDVI")
  );
}

function buildS2NdviCollection(start, end) {
  var col = ee.ImageCollection(dataSource)
    .filterBounds(aoi)
    .filterDate(start, end)
    .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", cloudMax))
    .map(maskS2)
    .map(scaleS2)
    .map(addMonth);

  if (useSeasonFilter) {
    col = col.filter(ee.Filter.inList("month", seasonMonths));
  }

  return col.map(addNdvi);
}

// ================= BUILD NDVI COMPOSITE =================
var s2Col = buildS2NdviCollection(recentStart, recentEnd);
print("S2 image count (filtered)", s2Col.size());

var ndviComposite = s2Col.select("NDVI").median().clip(aoi).rename("NDVI");
var obsCount = s2Col.select("NDVI").count().clip(aoi).rename("OBS_COUNT");

// Mask NDVI where there are too few observations
var ndvi = ndviComposite.updateMask(obsCount.gte(minObsCount));

// Optional analysis mask: remove built-up and water pixels so NDVI trend is less
// dominated by urban or permanent-water expansion around the mine.
var analysisMask = ee.Image(1).clip(aoi);
if (excludeBuiltUpAndWater) {
  var worldCover = ee.Image("ESA/WorldCover/v200/2021").clip(aoi);
  var builtUpMask = worldCover.eq(50); // built-up class
  var worldWaterMask = worldCover.eq(80); // permanent water class
  var jrcWaterMask = ee.Image("JRC/GSW1_4/GlobalSurfaceWater").select("occurrence").gte(50).clip(aoi);
  analysisMask = builtUpMask.not().and(worldWaterMask.not()).and(jrcWaterMask.not());
}

var ndviForStats = ndvi.updateMask(analysisMask);
var lowNdviMask = ndviForStats.lt(lowNdviThreshold).rename("LOW_NDVI");

// ================= BUFFER / RING GEOMETRIES =================
var distancesList = ee.List(bufferDistancesMeters).sort();

var cumulativeFeatures = ee.FeatureCollection(distancesList.map(function(d) {
  d = ee.Number(d);
  var geom = sitePoint.buffer(d);
  return ee.Feature(geom, {
    zone_type: "cumulative_buffer",
    radius_m: d,
    inner_m: 0,
    outer_m: d,
    label: ee.String("0-").cat(d.format()).cat(" m")
  });
}));

var ringFeatures = ee.FeatureCollection(distancesList.map(function(d) {
  d = ee.Number(d);
  var idx = distancesList.indexOf(d);
  var prev = ee.Number(ee.Algorithms.If(idx.eq(0), 0, distancesList.get(idx.subtract(1))));
  var outer = sitePoint.buffer(d);
  var inner = sitePoint.buffer(prev);
  var ring = outer.difference(inner, 1);
  return ee.Feature(ring, {
    zone_type: "ring",
    radius_m: d,
    inner_m: prev,
    outer_m: d,
    label: prev.format().cat("-").cat(d.format()).cat(" m")
  });
}));

Map.addLayer(cumulativeFeatures.style({
  color: "FFFFFF",
  fillColor: "00000000",
  width: 1
}), {}, "Cumulative buffers", false);

Map.addLayer(ringFeatures.style({
  color: "FFA500",
  fillColor: "00000000",
  width: 1
}), {}, "Rings", true);

// ================= ZONAL NDVI SUMMARY =================
function summarizeNdviByZones(zones, zoneName) {
  var stats = zones.map(function(f) {
    var obsQualifiedCount = ndvi.reduceRegion({
      reducer: ee.Reducer.count(),
      geometry: f.geometry(),
      scale: scaleMeters,
      bestEffort: true,
      maxPixels: 1e9
    }).get("NDVI");

    var analysisMaskCount = analysisMask.reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: f.geometry(),
      scale: scaleMeters,
      bestEffort: true,
      maxPixels: 1e9
    }).get("constant");

    var meanVal = ndviForStats.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: f.geometry(),
      scale: scaleMeters,
      bestEffort: true,
      maxPixels: 1e9
    }).get("NDVI");

    var countVal = ndviForStats.reduceRegion({
      reducer: ee.Reducer.count(),
      geometry: f.geometry(),
      scale: scaleMeters,
      bestEffort: true,
      maxPixels: 1e9
    }).get("NDVI");

    var validCount = ee.Number(ee.Algorithms.If(ee.Algorithms.IsEqual(countVal, null), 0, countVal));
    var obsQualified = ee.Number(ee.Algorithms.If(ee.Algorithms.IsEqual(obsQualifiedCount, null), 0, obsQualifiedCount));
    var analysisAvailable = ee.Number(ee.Algorithms.If(ee.Algorithms.IsEqual(analysisMaskCount, null), 0, analysisMaskCount));
    var zoneAreaM2 = ee.Number(f.geometry().area(1));
    var expectedPixels = zoneAreaM2.divide(scaleMeters * scaleMeters);
    var validPixelFraction = ee.Algorithms.If(
      expectedPixels.gt(0),
      validCount.divide(expectedPixels),
      null
    );

    var nullReason = ee.Algorithms.If(
      validCount.gt(0),
      "",
      ee.Algorithms.If(
        obsQualified.eq(0),
        "No pixels after minObsCount/time/cloud filtering",
        ee.Algorithms.If(
          ee.Algorithms.IsEqual(excludeBuiltUpAndWater, true),
          ee.Algorithms.If(
            analysisAvailable.eq(0),
            "All pixels removed by built-up/water analysis mask",
            "No pixels after combined masking"
          ),
          "No pixels after combined masking"
        )
      )
    );

    var lowNdviFraction = lowNdviMask.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: f.geometry(),
      scale: scaleMeters,
      bestEffort: true,
      maxPixels: 1e9
    }).get("LOW_NDVI");

    return f.set({
      mean: meanVal,
      count: countVal,
      site_name: siteName,
      metric: metricName,
      data_source: dataSource,
      composite_method: compositeMethod,
      analysis_window_start: recentStart.format("YYYY-MM-dd"),
      analysis_window_end: recentEnd.format("YYYY-MM-dd"),
      season_filter: useSeasonFilter,
      season_months: ee.String(useSeasonFilter ? seasonMonths.join(",") : ""),
      cloud_max: cloudMax,
      scale_m: scaleMeters,
      min_obs_count: minObsCount,
      obs_qualified_count: obsQualifiedCount,
      analysis_mask_count: analysisMaskCount,
      valid_pixel_fraction: validPixelFraction,
      exclude_built_up_and_water: excludeBuiltUpAndWater,
      low_ndvi_threshold: lowNdviThreshold,
      low_ndvi_fraction: lowNdviFraction,
      null_reason: nullReason
    });
  });

  print(zoneName + " NDVI zonal stats", stats);
  print(zoneName + " zones with null mean", stats.filter(ee.Filter.eq("mean", null)));
  return stats;
}

var cumulativeStats = summarizeNdviByZones(cumulativeFeatures, "Cumulative buffer");
var ringStats = summarizeNdviByZones(ringFeatures, "Ring");

// ================= VISUALIZATION =================
var visNdvi = {
  min: 0,
  max: 0.8,
  palette: ["8b4513", "f4d35e", "2ca25f", "006d2c"]
};
var visObs = {
  min: 0,
  max: 30,
  palette: ["2b2b2b", "f7f7f7", "00ff00"]
};

Map.addLayer(ndvi, visNdvi, "NDVI composite (recent)", true);
Map.addLayer(ndviForStats, visNdvi, "NDVI used for stats (optional land mask)", false);
Map.addLayer(obsCount, visObs, "QA: observation count", false);
Map.addLayer(analysisMask.selfMask(), {palette: ["C8E6C9"]}, "Analysis mask (built-up/water removed)", false);

// ================= CHARTS =================
var cumulativeStatsForChart = cumulativeStats.filter(ee.Filter.notNull(["mean"]));
var ringStatsForChart = ringStats.filter(ee.Filter.notNull(["mean"]));

print("Cumulative zones used in chart", cumulativeStatsForChart.size());
print("Ring zones used in chart", ringStatsForChart.size());

var cumulativeMeanChart = ui.Chart.feature.byFeature(cumulativeStatsForChart, "outer_m", ["mean"])
  .setChartType("LineChart")
  .setOptions({
    title: "NDVI Mean by Cumulative Buffer Distance",
    hAxis: {title: "Outer radius (m)"},
    vAxis: {title: "NDVI mean"},
    lineWidth: 2,
    pointSize: 4,
    colors: ["#2E7D32"]
  });
print(cumulativeMeanChart);

var ringMeanChart = ui.Chart.feature.byFeature(ringStatsForChart, "outer_m", ["mean"])
  .setChartType("ColumnChart")
  .setOptions({
    title: "NDVI Mean by Ring",
    hAxis: {title: "Ring outer radius (m)"},
    vAxis: {title: "NDVI mean"},
    colors: ["#66BB6A"]
  });
print(ringMeanChart);

// ================= OPTIONAL EXPORTS (UNCOMMENT IF NEEDED) =================
// Export.table.toDrive({
//   collection: cumulativeStats,
//   description: "san_sebastian_ndvi_cumulative_buffer_stats",
//   folder: "GEE_exports",
//   fileNamePrefix: "san_sebastian_ndvi_cumulative_buffer_stats",
//   fileFormat: "CSV"
// });
//
// Export.table.toDrive({
//   collection: ringStats,
//   description: "san_sebastian_ndvi_ring_stats",
//   folder: "GEE_exports",
//   fileNamePrefix: "san_sebastian_ndvi_ring_stats",
//   fileFormat: "CSV"
// });

print("NDVI progressive-buffer script ready. Change site coordinates, dates, and bufferDistancesMeters to reuse at other sites.");
