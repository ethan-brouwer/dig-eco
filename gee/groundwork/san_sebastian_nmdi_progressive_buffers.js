/*
  FILE: gee/groundwork/san_sebastian_nmdi_progressive_buffers.js
  PURPOSE: Lightweight visualization-first moisture screening around San Sebastian
           using progressive buffers/rings plus simple river context layers.

  DESIGN GOAL
  - Keep the workflow simple and easy to debug in the GEE Code Editor.
  - Visualize land moisture around the mine with NMDI.
  - Keep water separate from land moisture interpretation.
  - Show a river/water context layer without forcing downstream analysis yet.

  IMPORTANT LIMITATIONS
  - NMDI is more useful for surrounding land moisture/stress than for river water.
  - River visualization here is a screening/context layer, not a routed flow path.
  - This is a first-pass visualization script, not a full analysis workflow.
*/

// ================= USER SETTINGS =================
var siteName = "San Sebastian Mine (MRDS)";
var siteLon = -87.93002;
var siteLat = 13.6509;

// Progressive distances from site center (meters).
var bufferDistancesMeters = [100, 250, 500, 1000, 1500];

// Temporal / QA settings
var startDate = "2025-11-01";
var endDate = "2026-03-01";   // exclusive end date
var cloudMax = 60;
var seasonMonths = [11, 12, 1, 2, 3, 4];
var useSeasonFilter = true;

// Visualization / screening settings
var scaleMeters = 20;
var minObsCount = 2;
var jrcWaterOccurrenceThreshold = 50;
var streamAccThreshold = 20;

// ================= MAP SETUP =================
Map.setOptions("SATELLITE");

var sitePoint = ee.Geometry.Point([siteLon, siteLat]);
var maxRadius = ee.Number(ee.List(bufferDistancesMeters).reduce(ee.Reducer.max()));
var aoi = sitePoint.buffer(maxRadius);

Map.centerObject(sitePoint, 13);
Map.addLayer(sitePoint, {color: "FF0000"}, siteName, true);
Map.addLayer(aoi, {color: "FFFFFF"}, "AOI (max buffer)", false);

// ================= DATE WINDOW =================
var recentStart = ee.Date(startDate);
var recentEnd = ee.Date(endDate);

print("Analysis window", recentStart.format("YYYY-MM-dd"), recentEnd.format("YYYY-MM-dd"));
print("Buffer distances (m)", bufferDistancesMeters);
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
  var optical = img.select(["B2", "B3", "B4", "B8", "B11", "B12"]).multiply(0.0001);
  return img.addBands(optical, null, true);
}

function addMonth(img) {
  return img.set("month", ee.Date(img.get("system:time_start")).get("month"));
}

function addIndices(img) {
  var green = img.select("B3");
  var red = img.select("B4");
  var nir = img.select("B8");
  var swir1 = img.select("B11");
  var swir2 = img.select("B12");

  var ndwi = img.normalizedDifference(["B3", "B8"]).rename("NDWI");
  var mndwi = img.normalizedDifference(["B3", "B11"]).rename("MNDWI");
  var ndti = img.normalizedDifference(["B4", "B3"]).rename("NDTI");
  var tssProxy = red.multiply(1000).rename("TSS_PROXY");
  var redGreen = red.divide(green).rename("RED_GREEN");
  var nmdi = img.expression(
    "(NIR - (SWIR1 - SWIR2)) / (NIR + (SWIR1 - SWIR2))",
    {NIR: nir, SWIR1: swir1, SWIR2: swir2}
  ).rename("NMDI");

  return img.addBands([ndwi, mndwi, ndti, tssProxy, redGreen, nmdi]);
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

  return col.map(addIndices);
}

// ================= BUILD COMPOSITE =================
var s2Col = buildS2Collection(recentStart, recentEnd);
print("S2 image count (filtered)", s2Col.size());

var recent = s2Col.median().clip(aoi);
var obsCount = s2Col.select("B8").count().clip(aoi).rename("OBS_COUNT");
var recentTrueColor = recent.select(["B4", "B3", "B2"]).rename(["R", "G", "B"]);
var recentFalseColor = recent.select(["B8", "B4", "B3"]).rename(["NIR", "R", "G"]);

// ================= WATER / RIVER CONTEXT =================
var jrcWater = ee.Image("JRC/GSW1_4/GlobalSurfaceWater")
  .select("occurrence")
  .unmask(0)
  .clip(aoi);
var permanentWater = jrcWater.gte(jrcWaterOccurrenceThreshold).rename("PERMANENT_WATER");

function waterMask(img) {
  var opticalWater = img.select("NDWI").gt(0.12).or(img.select("MNDWI").gt(0.15));
  return opticalWater;
}

var currentWater = waterMask(recent).selfMask().rename("CURRENT_WATER");
var hydroAcc = ee.Image("WWF/HydroSHEDS/15ACC").clip(aoi);
var streamMask = hydroAcc.gte(streamAccThreshold)
  .focal_max(30, "circle", "meters")
  .selfMask()
  .rename("STREAM_MASK");
var riverContext = permanentWater.selfMask().or(streamMask).selfMask().rename("RIVER_CONTEXT");
var riverWaterMask = currentWater.and(riverContext).selfMask().rename("RIVER_WATER_MASK");

// Land-focused NMDI view: keep river/water separate from land moisture.
var nmdi = recent.select("NMDI")
  .updateMask(obsCount.gte(minObsCount))
  .updateMask(riverWaterMask.not())
  .rename("NMDI");

// Water-focused river condition views.
var ndtiWater = recent.select("NDTI")
  .updateMask(obsCount.gte(minObsCount))
  .updateMask(riverWaterMask)
  .rename("NDTI_WATER");

var tssWater = recent.select("TSS_PROXY")
  .updateMask(obsCount.gte(minObsCount))
  .updateMask(riverWaterMask)
  .rename("TSS_WATER");

var redGreenWater = recent.select("RED_GREEN")
  .updateMask(obsCount.gte(minObsCount))
  .updateMask(riverWaterMask)
  .rename("RED_GREEN_WATER");

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

print("River water pixel count", riverWaterMask.reduceRegion({
  reducer: ee.Reducer.count(),
  geometry: aoi,
  scale: scaleMeters,
  bestEffort: true,
  maxPixels: 1e9
}).values().get(0));

print("River NDTI range", ndtiWater.reduceRegion({
  reducer: ee.Reducer.percentile([5, 50, 95]),
  geometry: aoi,
  scale: scaleMeters,
  bestEffort: true,
  maxPixels: 1e9
}));

print("River TSS proxy range", tssWater.reduceRegion({
  reducer: ee.Reducer.percentile([5, 50, 95]),
  geometry: aoi,
  scale: scaleMeters,
  bestEffort: true,
  maxPixels: 1e9
}));

// ================= VISUALIZATION =================
var visNmdi = {min: -0.6, max: 0.6, palette: ["8c510a", "f6e8c3", "c7eae5", "01665e"]};
var visNdti = {min: -0.08, max: 0.08, palette: ["2166AC", "F7F7F7", "B2182B"]};
var visTss = {min: 40, max: 180, palette: ["08306B", "41B6C4", "FFFFBF", "FDAE61", "D73027"]};
var visRedGreen = {min: 0.7, max: 1.4, palette: ["F7FBFF", "FDAE61", "D73027"]};
var visObsCount = {min: 0, max: 20, palette: ["2b2b2b", "f7f7f7", "00ff00"]};

Map.addLayer(recentTrueColor, {bands: ["R", "G", "B"], min: 0.02, max: 0.30}, "Recent composite (true color RGB)", true);
Map.addLayer(recentFalseColor, {bands: ["NIR", "R", "G"], min: 0.03, max: 0.40}, "Recent composite (false color NIR-R-G)", false);
Map.addLayer(nmdi, visNmdi, "NMDI (land only; water masked)", true);
Map.addLayer(permanentWater.selfMask(), {palette: ["4FC3F7"]}, "River context: JRC permanent water", false);
Map.addLayer(streamMask, {palette: ["00FFFF"]}, "River context: HydroSHEDS stream mask", false);
Map.addLayer(riverContext, {palette: ["80DEEA"]}, "River context: combined water corridor", true);
Map.addLayer(currentWater, {palette: ["64B5F6"]}, "Current optical water mask", false);
Map.addLayer(riverWaterMask, {palette: ["00BFFF"]}, "River analysis mask (current water in corridor)", true);
Map.addLayer(ndtiWater, visNdti, "Sediment proxy: NDTI (water only)", true);
Map.addLayer(tssWater, visTss, "Sediment proxy: red reflectance / TSS proxy (water only)", false);
Map.addLayer(redGreenWater, visRedGreen, "Sediment proxy: red-green ratio (water only)", false);
Map.addLayer(obsCount, visObsCount, "QA observation count", false);

Map.addLayer(cumulativeFeatures.style({
  color: "FFFFFF",
  fillColor: "00000000",
  width: 1
}), {}, "Cumulative buffers", false);

Map.addLayer(ringFeatures.style({
  color: "FFA500",
  fillColor: "00000000",
  width: 2
}), {}, "Progressive rings", true);

print("NMDI is being used here for surrounding land moisture/stress, not for the river itself.");
print("River-condition layers are water-only screening proxies: NDTI, red reflectance/TSS proxy, and red-green ratio.");
print("This simplified version uses a shorter date range, smaller AOI, and a looser river-water mask so the river layers render more reliably.");
