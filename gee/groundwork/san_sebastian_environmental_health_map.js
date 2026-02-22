/*
  FILE: gee/groundwork/san_sebastian_environmental_health_map.js
  PURPOSE: San Sebastian Mine (El Salvador) "current conditions" map layers
           for environmental health screening in the GEE Code Editor.
  INPUTS: Sentinel-2 SR Harmonized, JRC Global Surface Water, optional ESA WorldCover.
  OUTPUTS: Map layers + checklist panel (no export by default).

  IMPORTANT LIMITATION
  - This script does NOT directly detect sulfide, nickel, gold, or dissolved metals.
  - It visualizes optical/surface proxies that can indicate mine disturbance,
    vegetation stress, sediment/turbidity, and iron-oxide-rich materials.
  - Use in-situ chemistry/lab data to confirm contamination.
*/

// ================= USER SETTINGS =================
var mineName = "San Sebastian Mine (MRDS)";
var mineLon = -87.93002;  // MRDS CSV in repo
var mineLat = 13.6509;

var aoiRadiusMeters = 12000;
var innerFocusMeters = 3000;
var recentMonths = 12;      // "current" window; keeps at least one dry-season pass most of the year
var baselineYears = 5;      // prior multi-year baseline for anomaly comparison
var cloudMax = 60;
var seasonMonths = [11, 12, 1, 2, 3, 4]; // drier season often better for turbidity visibility

// ================= MAP SETUP =================
Map.setOptions("SATELLITE");

var minePoint = ee.Geometry.Point([mineLon, mineLat]);
var aoi = minePoint.buffer(aoiRadiusMeters);
var innerFocus = minePoint.buffer(innerFocusMeters);

Map.centerObject(minePoint, 12);
Map.addLayer(aoi, {color: "FFFFFF"}, "AOI (12 km)", false);
Map.addLayer(innerFocus, {color: "FFAA00"}, "Inner focus (3 km)", false);
Map.addLayer(minePoint, {color: "FF0000"}, mineName, true);

// ================= DATE WINDOWS =================
var now = ee.Date(Date.now());
var recentStart = now.advance(-recentMonths, "month");
var recentEnd = now;

// Baseline excludes the most recent year to reduce overlap with "current".
var baselineStart = now.advance(-(baselineYears + 1), "year");
var baselineEnd = now.advance(-1, "year");

print("Recent window", recentStart.format("YYYY-MM-dd"), recentEnd.format("YYYY-MM-dd"));
print("Baseline window", baselineStart.format("YYYY-MM-dd"), baselineEnd.format("YYYY-MM-dd"));

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
  // Reflectance scale factor for S2 SR bands.
  var optical = img.select([
    "B2", "B3", "B4", "B8", "B11", "B12"
  ]).multiply(0.0001);
  return img.addBands(optical, null, true);
}

function addMonth(img) {
  var month = ee.Date(img.get("system:time_start")).get("month");
  return img.set("month", month);
}

function addIndices(img) {
  var blue = img.select("B2");
  var green = img.select("B3");
  var red = img.select("B4");
  var nir = img.select("B8");
  var swir1 = img.select("B11");
  var swir2 = img.select("B12");

  var ndvi = img.normalizedDifference(["B8", "B4"]).rename("NDVI");
  var savi = img.expression(
    "((NIR - RED) / (NIR + RED + L)) * (1 + L)",
    {NIR: nir, RED: red, L: 0.5}
  ).rename("SAVI");
  var ndmi = img.normalizedDifference(["B8", "B11"]).rename("NDMI");
  var ndwi = img.normalizedDifference(["B3", "B8"]).rename("NDWI");
  var mndwi = img.normalizedDifference(["B3", "B11"]).rename("MNDWI");
  var ndti = img.normalizedDifference(["B4", "B3"]).rename("NDTI"); // sediment / water discoloration proxy

  var bsi = img.expression(
    "((SWIR + RED) - (NIR + BLUE)) / ((SWIR + RED) + (NIR + BLUE))",
    {SWIR: swir1, RED: red, NIR: nir, BLUE: blue}
  ).rename("BSI");

  // Mining/mineral surface proxies commonly used in multispectral screening.
  var ioi = red.divide(blue).rename("IOI");           // iron oxide index proxy
  var ferrous = swir1.divide(nir).rename("FERROUS");  // ferrous iron proxy ratio
  var clay = swir1.divide(swir2).rename("CLAY");      // clay/alteration proxy ratio
  var kaolinite = swir1.divide(swir2).rename("KAOLINITE");

  // Water turbidity / suspended sediment proxies (optical only).
  var redGreen = red.divide(green).rename("RED_GREEN");
  var redNir = red.divide(nir).rename("RED_NIR");
  var tssProxy = red.multiply(1000).rename("TSS_PROXY"); // relative proxy only

  return img.addBands([
    ndvi, savi, ndmi, ndwi, mndwi, ndti, bsi,
    ioi, ferrous, clay, kaolinite,
    redGreen, redNir, tssProxy
  ]);
}

function buildS2Collection(start, end) {
  return ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
    .filterBounds(aoi)
    .filterDate(start, end)
    .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", cloudMax))
    .map(maskS2)
    .map(scaleS2)
    .map(addMonth)
    .filter(ee.Filter.inList("month", seasonMonths))
    .map(addIndices);
}

function compositeWithCounts(col, label) {
  print(label + " image count", col.size());
  return col.median().clip(aoi);
}

var recentCol = buildS2Collection(recentStart, recentEnd);
var baselineCol = buildS2Collection(baselineStart, baselineEnd);

var recent = compositeWithCounts(recentCol, "Recent");
var baseline = compositeWithCounts(baselineCol, "Baseline");

// ================= WATER MASK + TURBIDITY =================
var jrcWater = ee.Image("JRC/GSW1_4/GlobalSurfaceWater").select("occurrence");
var permanentWater = jrcWater.gte(50);

function waterMask(img) {
  var water = img.select("NDWI").gt(0.12).or(img.select("MNDWI").gt(0.15));
  return water.and(permanentWater);
}

var recentWaterMask = waterMask(recent);
var baselineWaterMask = waterMask(baseline);
var stableWaterMask = recentWaterMask.and(baselineWaterMask);

var recentTurbidity = recent.select("TSS_PROXY").updateMask(stableWaterMask).rename("TSS_RECENT");
var baselineTurbidity = baseline.select("TSS_PROXY").updateMask(stableWaterMask).rename("TSS_BASELINE");
var turbidityAnomaly = recentTurbidity.subtract(baselineTurbidity).rename("TSS_ANOMALY");
var recentNDTIWater = recent.select("NDTI").updateMask(stableWaterMask).rename("NDTI_WATER");

// ================= SCREENING SCORES (RED = HIGHER CONCERN) =================
function norm(img, min, max) {
  return img.unitScale(min, max).clamp(0, 1);
}

function invertNorm(img, min, max) {
  return ee.Image(1).subtract(norm(img, min, max));
}

var landMask = stableWaterMask.not();
var waterMaskImg = stableWaterMask;

// Land screening: exposed/disturbed surfaces + low vegetation + iron/mineral proxies
var landStress = ee.ImageCollection.fromImages([
  norm(recent.select("BSI"), 0.0, 0.45),
  invertNorm(recent.select("NDVI"), 0.2, 0.8),
  invertNorm(recent.select("NDMI"), -0.1, 0.4),
  norm(recent.select("IOI"), 0.8, 2.0),
  norm(recent.select("FERROUS"), 0.8, 1.8),
  norm(recent.select("CLAY"), 0.8, 1.6)
]).mean().updateMask(landMask).rename("LAND_STRESS");

// Water screening: recent turbidity + change from baseline + NDTI water discoloration
var waterStress = ee.ImageCollection.fromImages([
  norm(recentTurbidity, 20, 250),
  norm(turbidityAnomaly, 0, 80),   // only positive anomalies emphasized
  norm(recentNDTIWater, 0.0, 0.3)
]).mean().updateMask(waterMaskImg).rename("WATER_STRESS");

var combinedScreen = landStress.unmask(0)
  .max(waterStress.unmask(0))
  .rename("MINE_IMPACT_SCREEN");

var hotSpots = combinedScreen.gt(0.65).selfMask();

// ================= VISUALIZATION PRESETS =================
var visTrueColor = {bands: ["B4", "B3", "B2"], min: 0.02, max: 0.30};
var visFalseColor = {bands: ["B8", "B4", "B3"], min: 0.03, max: 0.45};

var visNdvi = {min: 0, max: 0.8, palette: ["8b4513", "f4d35e", "2ca25f", "006d2c"]};
var visBsi = {min: -0.3, max: 0.5, palette: ["08306b", "f7f7f7", "ff0000"]};   // red = more bare/disturbed
var visNdmi = {min: -0.4, max: 0.5, palette: ["a6611a", "f7f7f7", "1a9850"]};
var visNDTI = {min: -0.3, max: 0.4, palette: ["2166ac", "f7f7f7", "b2182b"]};  // red = more reddish/sediment signal
var visRatio = {min: 0.7, max: 2.0, palette: ["f7fbff", "fdae61", "d73027"]};
var visTurbidity = {min: 20, max: 300, palette: ["08306b", "41b6c4", "ffffbf", "fdae61", "d73027"]};
var visTurbidityAnom = {min: -80, max: 80, palette: ["2166ac", "f7f7f7", "b2182b"]};
var visStress = {min: 0, max: 1, palette: ["ffffcc", "fd8d3c", "bd0026"]};      // red = higher screening concern

// ================= MAP LAYERS =================
Map.addLayer(recent, visTrueColor, "Recent composite (true color)", true);
Map.addLayer(recent, visFalseColor, "Recent composite (false color NIR)", false);

// Core health layers requested
Map.addLayer(recent.select("NDVI"), visNdvi, "NDVI (vegetation health)", true);
Map.addLayer(recent.select("SAVI"), visNdvi, "SAVI (soil-adjusted veg)", false);
Map.addLayer(recent.select("BSI"), visBsi, "BSI (bare/disturbed soil, red=high)", true);

// Additional mine-impact proxies (literature-style multispectral screening)
Map.addLayer(recent.select("NDMI"), visNdmi, "NDMI (moisture stress)", false);
Map.addLayer(recent.select("NDTI"), visNDTI, "NDTI (sediment/discoloration proxy)", false);
Map.addLayer(recent.select("IOI"), visRatio, "Iron Oxide Index proxy (red=high)", true);
Map.addLayer(recent.select("FERROUS"), {min: 0.7, max: 1.8, palette: ["ffffff", "fdae61", "d73027"]}, "Ferrous proxy ratio (red=high)", false);
Map.addLayer(recent.select("CLAY"), {min: 0.8, max: 1.6, palette: ["f7fcf0", "addd8e", "006837"]}, "Clay/alteration proxy", false);

// Water and turbidity
Map.addLayer(stableWaterMask.selfMask(), {palette: ["00BFFF"]}, "Stable water mask", false);
Map.addLayer(recentTurbidity, visTurbidity, "Turbidity/TSS proxy (water)", true);
Map.addLayer(turbidityAnomaly, visTurbidityAnom, "Turbidity anomaly vs baseline (red=increased)", true);
Map.addLayer(recentNDTIWater, {min: -0.1, max: 0.3, palette: ["2166ac", "f7f7f7", "d73027"]}, "Water NDTI (discoloration)", false);

// Screening layers (red-forward for quick issue scan)
Map.addLayer(landStress, visStress, "Land mine-impact screening (red=high)", true);
Map.addLayer(waterStress, visStress, "Water contamination screening (red=high)", true);
Map.addLayer(combinedScreen, visStress, "Combined environmental screening (red=high)", true);
Map.addLayer(hotSpots, {palette: ["FF0000"]}, "Hotspots > 0.65", true);

// Optional contextual layer: current bare ground cover class
var worldCover = ee.Image("ESA/WorldCover/v200/2021");
var bareCover = worldCover.eq(60).clip(aoi).selfMask();
Map.addLayer(bareCover, {palette: ["FFD54F"]}, "WorldCover bare/sparse (2021)", false);

// ================= QUICK CHECK PANEL =================
var panel = ui.Panel({
  style: {position: "top-right", width: "360px", padding: "8px"}
});

panel.add(ui.Label({
  value: "San Sebastian Environmental Health Layers",
  style: {fontWeight: "bold", fontSize: "13px"}
}));

panel.add(ui.Label("Direct-ish optical layers (current surface condition):"));
panel.add(ui.Label("[x] NDVI / SAVI vegetation condition"));
panel.add(ui.Label("[x] BSI exposed/bare soil disturbance"));
panel.add(ui.Label("[x] NDMI moisture stress"));
panel.add(ui.Label("[x] Turbidity/TSS proxy over water"));

panel.add(ui.Label("Mining-related screening proxies (not direct chemistry):"));
panel.add(ui.Label("[x] NDTI water/soil discoloration"));
panel.add(ui.Label("[x] Iron oxide index proxy (IOI)"));
panel.add(ui.Label("[x] Ferrous proxy ratio"));
panel.add(ui.Label("[x] Clay/alteration proxy ratio"));
panel.add(ui.Label("[x] Baseline-vs-current turbidity anomaly"));
panel.add(ui.Label("[x] Red hotspot screening layers"));

panel.add(ui.Label({
  value: "Red means higher screening concern (disturbance/stress/turbidity proxy), not confirmed sulfide/nickel concentration.",
  style: {fontSize: "11px", color: "B22222"}
}));

Map.add(panel);

print("NOTE: Sulfide, nickel, gold, and dissolved metals require water/soil chemistry sampling.");
print("This script is for remote sensing triage and hotspot targeting.");

// ================= OPTIONAL EXPORTS (UNCOMMENT IF NEEDED) =================
// Export.image.toDrive({
//   image: recent.select(["NDVI", "BSI", "NDMI", "IOI", "FERROUS", "CLAY"]),
//   description: "san_sebastian_recent_land_indices",
//   folder: "GEE_exports",
//   fileNamePrefix: "san_sebastian_recent_land_indices",
//   region: aoi,
//   scale: 10,
//   maxPixels: 1e13
// });
//
// Export.image.toDrive({
//   image: ee.Image.cat([recentTurbidity, baselineTurbidity, turbidityAnomaly, waterStress]),
//   description: "san_sebastian_turbidity_screening",
//   folder: "GEE_exports",
//   fileNamePrefix: "san_sebastian_turbidity_screening",
//   region: aoi,
//   scale: 10,
//   maxPixels: 1e13
// });
