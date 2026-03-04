/*
  FILE: gee/groundwork/san_sebastian_site_characterization_simple.js
  PURPOSE: Lightweight site characterization for San Sebastian Mine (El Salvador)
           using recent Sentinel-2 optical proxies in the GEE Code Editor.

  NOTES
  - This is a simplified "current conditions" characterization script.
  - It omits baseline anomaly scoring, unified risk/health scores, river corridor
    splitting, and monthly time-series charts from the advanced script.
  - All layers are screening proxies only, not direct contamination measurements.
*/

// ================= USER SETTINGS =================
var mineName = "San Sebastian Mine (MRDS)";
var mineLon = -87.93002;
var mineLat = 13.6509;

var aoiRadiusMeters = 3000;
var recentMonths = 12;
var cloudMax = 60;
var seasonMonths = [11, 12, 1, 2, 3, 4];
var useSeasonFilter = true;

// ================= MAP SETUP =================
Map.setOptions("SATELLITE");

var minePoint = ee.Geometry.Point([mineLon, mineLat]);
var aoi = minePoint.buffer(aoiRadiusMeters);

Map.centerObject(minePoint, 12);
Map.addLayer(aoi, {color: "FFFFFF"}, "AOI (3 km)", true);
Map.addLayer(minePoint, {color: "FF0000"}, mineName, true);

// ================= DATE WINDOW =================
var now = ee.Date(Date.now());
var recentStart = now.advance(-recentMonths, "month");
var recentEnd = now;
print("Recent window", recentStart.format("YYYY-MM-dd"), recentEnd.format("YYYY-MM-dd"));

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
  var blue = img.select("B2");
  var green = img.select("B3");
  var red = img.select("B4");
  var nir = img.select("B8");
  var swir1 = img.select("B11");

  var ndvi = img.normalizedDifference(["B8", "B4"]).rename("NDVI");
  var ndmi = img.normalizedDifference(["B8", "B11"]).rename("NDMI");
  var ndwi = img.normalizedDifference(["B3", "B8"]).rename("NDWI");
  var mndwi = img.normalizedDifference(["B3", "B11"]).rename("MNDWI");
  var ndti = img.normalizedDifference(["B4", "B3"]).rename("NDTI");
  var redGreen = red.divide(green).rename("RED_GREEN");

  var bsi = img.expression(
    "((SWIR + RED) - (NIR + BLUE)) / ((SWIR + RED) + (NIR + BLUE))",
    {SWIR: swir1, RED: red, NIR: nir, BLUE: blue}
  ).rename("BSI");

  var ioi = red.divide(blue).rename("IOI");
  var tssProxy = red.multiply(1000).rename("TSS_PROXY");

  return img.addBands([ndvi, ndmi, ndwi, mndwi, ndti, redGreen, bsi, ioi, tssProxy]);
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

var recentCol = buildS2Collection(recentStart, recentEnd);
print("Recent image count", recentCol.size());

var recent = recentCol.median().clip(aoi);
var recentObsCount = recentCol.select("B4").count().clip(aoi).rename("RECENT_OBS_COUNT");

// ================= RIVER CONTEXT + WATER PROXIES =================
var jrcWater = ee.Image("JRC/GSW1_4/GlobalSurfaceWater").select("occurrence");
var permanentWater = jrcWater.gte(50).clip(aoi).rename("PERMANENT_WATER");
var hydroAcc = ee.Image("WWF/HydroSHEDS/15ACC").clip(aoi).rename("HYDRO_ACC");
var coarseRiverContext = hydroAcc.gte(50).selfMask().rename("COARSE_RIVER_CONTEXT");

function waterMask(img) {
  var opticalWater = img.select("NDWI").gt(0.12).or(img.select("MNDWI").gt(0.15));
  return opticalWater.and(permanentWater);
}

var currentWaterMask = waterMask(recent).selfMask().rename("CURRENT_WATER");
var riverAnalysisMask = currentWaterMask.and(coarseRiverContext.unmask(0)).selfMask().rename("RIVER_ANALYSIS_MASK");
var excludedFromRiverMask = currentWaterMask.and(riverAnalysisMask.not()).selfMask().rename("EXCLUDED_FROM_RIVER_MASK");

var recentTurbidity = recent.select("TSS_PROXY").updateMask(riverAnalysisMask).rename("TSS_RECENT");
var recentNDTIWater = recent.select("NDTI").updateMask(riverAnalysisMask).rename("NDTI_WATER");
var recentRedGreenWater = recent.select("RED_GREEN").updateMask(riverAnalysisMask).rename("RED_GREEN_WATER");

var waterProxyComposite = ee.Image.cat([
  recent.select("TSS_PROXY").unitScale(20, 220).clamp(0, 1),
  recent.select("NDTI").unitScale(-0.1, 0.25).clamp(0, 1),
  recent.select("RED_GREEN").unitScale(0.7, 2.0).clamp(0, 1)
]).updateMask(riverAnalysisMask).rename(["TSS_RGB", "NDTI_RGB", "RED_GREEN_RGB"]);

// ================= VISUALIZATION =================
var visTrueColor = {bands: ["B4", "B3", "B2"], min: 0.02, max: 0.30};
var visFalseColor = {bands: ["B8", "B4", "B3"], min: 0.03, max: 0.45};
var visNdvi = {min: 0, max: 0.8, palette: ["8b4513", "f4d35e", "2ca25f", "006d2c"]};
var visNdmi = {min: -0.4, max: 0.5, palette: ["a6611a", "f7f7f7", "1a9850"]};
var visBsi = {min: -0.3, max: 0.5, palette: ["08306b", "f7f7f7", "ff0000"]};
var visNDTI = {min: -0.3, max: 0.4, palette: ["2166ac", "f7f7f7", "b2182b"]};
var visRatio = {min: 0.7, max: 2.0, palette: ["f7fbff", "fdae61", "d73027"]};
var visTurbidity = {min: 20, max: 300, palette: ["08306b", "41b6c4", "ffffbf", "fdae61", "d73027"]};
var visObsCount = {min: 0, max: 30, palette: ["2b2b2b", "f7f7f7", "00ff00"]};
var visProxyComposite = {bands: ["TSS_RGB", "NDTI_RGB", "RED_GREEN_RGB"], min: 0, max: 1};

// ================= MAP LAYERS =================
Map.addLayer(recent, visTrueColor, "Recent composite (true color)", true);
Map.addLayer(recent, visFalseColor, "Recent composite (false color NIR)", false);

Map.addLayer(recent.select("NDVI"), visNdvi, "NDVI (vegetation condition)", true);
Map.addLayer(recent.select("NDMI"), visNdmi, "NDMI (moisture stress proxy)", false);
Map.addLayer(recent.select("BSI"), visBsi, "BSI (bare/disturbed surface proxy)", true);
Map.addLayer(recent.select("IOI"), visRatio, "IOI (iron-oxide proxy)", false);

Map.addLayer(coarseRiverContext, {palette: ["90CAF9"]}, "Coarse river context (HydroSHEDS)", false);
Map.addLayer(permanentWater.selfMask(), {palette: ["4FC3F7"]}, "Permanent water (JRC)", false);
Map.addLayer(currentWaterMask, {palette: ["00BFFF"]}, "Current water mask (optical + JRC)", false);
Map.addLayer(riverAnalysisMask, {palette: ["FFD54F"]}, "River analysis mask", true);
Map.addLayer(excludedFromRiverMask, {palette: ["FF00FF"]}, "QA excluded from river analysis mask", false);
Map.addLayer(waterProxyComposite, visProxyComposite, "Water proxy composite (TSS / NDTI / red-green)", true);
Map.addLayer(recentTurbidity, visTurbidity, "Water turbidity/TSS proxy", true);
Map.addLayer(recentNDTIWater, visNDTI, "Water NDTI (discoloration proxy)", false);
Map.addLayer(recentRedGreenWater, visRatio, "Water red/green ratio", false);

Map.addLayer(recentObsCount, visObsCount, "QA: recent observation count", false);

// Optional context layer
var worldCover = ee.Image("ESA/WorldCover/v200/2021");
var bareCover = worldCover.eq(60).clip(aoi).selfMask();
Map.addLayer(bareCover, {palette: ["FFD54F"]}, "WorldCover bare/sparse (2021)", false);

// ================= SIMPLE PANEL =================
var panel = ui.Panel({
  style: {position: "top-right", width: "330px", padding: "8px"}
});

panel.add(ui.Label({
  value: "San Sebastian Site Characterization (Simple)",
  style: {fontWeight: "bold", fontSize: "13px"}
}));
panel.add(ui.Label("Current-condition optical screening layers:"));
panel.add(ui.Label("[x] Recent Sentinel-2 composite"));
panel.add(ui.Label("[x] NDVI (vegetation condition)"));
panel.add(ui.Label("[x] BSI (bare/disturbed surface proxy)"));
panel.add(ui.Label("[x] NDMI (moisture stress proxy)"));
panel.add(ui.Label("[x] IOI (iron-oxide proxy)"));
panel.add(ui.Label("[x] Coarse river context + river analysis mask"));
panel.add(ui.Label("[x] Water proxy composite (TSS / NDTI / red-green)"));
panel.add(ui.Label("[x] Water mask + turbidity/NDTI water proxies"));
panel.add(ui.Label("[x] QA image-count layer"));
panel.add(ui.Label({
  value: "Proxies only: use field/lab sampling to confirm metals or contamination.",
  style: {fontSize: "11px", color: "B22222"}
}));

Map.add(panel);

print("NOTE: Simplified site characterization mode (no baseline anomaly scoring or unified risk/health score).");
print("Use the advanced script for anomaly heatmaps, significance masks, and upstream/downstream analysis.");
