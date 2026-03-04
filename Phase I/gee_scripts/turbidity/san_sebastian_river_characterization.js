/*
  FILE: Phase I/gee_scripts/turbidity/san_sebastian_river_characterization.js
  PURPOSE: Visualization-first river characterization for San Sebastian Mine.

  DESIGN GOAL
  - Focus only on identifying the river analysis area and water-only proxies.
  - Keep the masking logic transparent and easy to debug in the GEE Code Editor.
  - Avoid inferred upstream/downstream bins at this stage.

  IMPORTANT LIMITATIONS
  - This is a screening script, not a contamination measurement workflow.
  - JRC is coarse and can over-widen the channel footprint.
  - Optical water masks can miss narrow river segments under cloud/shadow,
    riparian cover, shallow water, or low-contrast conditions.
*/

// ================= USER SETTINGS =================
var siteName = "San Sebastian Mine (MRDS)";
var siteLon = -87.93002;
var siteLat = 13.6509;

var aoiRadiusMeters = 5000;
var recentMonths = 12;
var cloudMax = 60;
var useSeasonFilter = false;
var seasonMonths = [11, 12, 1, 2, 3, 4];

var scaleMeters = 20;
var minObsCount = 2;

// River-mask controls.
var jrcOccurrenceThreshold = 5;
var jrcCorridorExpandMeters = 60;
var hydroAccThreshold = 20;
var hydroCorridorExpandMeters = 60;
var ndwiThreshold = 0.05;
var mndwiThreshold = 0.02;

// ================= MAP SETUP =================
Map.setOptions("SATELLITE");

var sitePoint = ee.Geometry.Point([siteLon, siteLat]);
var aoi = sitePoint.buffer(aoiRadiusMeters);

Map.centerObject(sitePoint, 13);
Map.addLayer(sitePoint, {color: "FF0000"}, siteName, true);
Map.addLayer(aoi, {color: "FFFFFF"}, "AOI (5 km)", false);

// ================= DATE WINDOW =================
var recentEnd = ee.Date(Date.now());
var recentStart = recentEnd.advance(-recentMonths, "month");

print("Recent window", recentStart.format("YYYY-MM-dd"), recentEnd.format("YYYY-MM-dd"));
print("JRC occurrence threshold", jrcOccurrenceThreshold);
print("JRC corridor expansion (m)", jrcCorridorExpandMeters);
print("HydroSHEDS accumulation threshold", hydroAccThreshold);
print("HydroSHEDS corridor expansion (m)", hydroCorridorExpandMeters);
print("NDWI threshold", ndwiThreshold);
print("MNDWI threshold", mndwiThreshold);

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

// ================= COMPOSITE =================
var recentCol = buildS2Collection(recentStart, recentEnd);
var recent = recentCol.median().clip(aoi);
var obsCount = recentCol.select("B4").count().clip(aoi).rename("OBS_COUNT");

print("Recent image count", recentCol.size());

// ================= RIVER CHARACTERIZATION =================
var jrcOccurrence = ee.Image("JRC/GSW1_4/GlobalSurfaceWater")
  .select("occurrence")
  .unmask(0)
  .clip(aoi)
  .rename("JRC_OCCURRENCE");

var jrcRiverSeed = jrcOccurrence.gte(jrcOccurrenceThreshold)
  .selfMask()
  .rename("JRC_RIVER_SEED");

var jrcRiverCorridor = jrcRiverSeed
  .focal_max(jrcCorridorExpandMeters, "circle", "meters")
  .selfMask()
  .rename("JRC_RIVER_CORRIDOR");

var hydroAcc = ee.Image("WWF/HydroSHEDS/15ACC")
  .clip(aoi)
  .rename("HYDRO_ACC");

var hydroRiverSeed = hydroAcc.gte(hydroAccThreshold)
  .selfMask()
  .rename("HYDRO_RIVER_SEED");

var hydroRiverCorridor = hydroRiverSeed
  .focal_max(hydroCorridorExpandMeters, "circle", "meters")
  .selfMask()
  .rename("HYDRO_RIVER_CORRIDOR");

var opticalWater = recent.select("NDWI").gt(ndwiThreshold)
  .or(recent.select("MNDWI").gt(mndwiThreshold))
  .updateMask(obsCount.gte(minObsCount))
  .selfMask()
  .rename("OPTICAL_WATER");

var riverAnalysisMask = opticalWater
  .selfMask()
  .rename("RIVER_ANALYSIS_MASK");

var hydroConstrainedRiverMask = opticalWater
  .and(hydroRiverCorridor.unmask(0))
  .selfMask()
  .rename("HYDRO_CONSTRAINED_RIVER_MASK");

var jrcConstrainedRiverMask = opticalWater
  .and(jrcRiverCorridor.unmask(0))
  .selfMask()
  .rename("JRC_CONSTRAINED_RIVER_MASK");

var opticalOutsideHydro = opticalWater
  .and(hydroRiverCorridor.unmask(0).not())
  .selfMask()
  .rename("OPTICAL_OUTSIDE_HYDRO");

var opticalOutsideJrc = opticalWater
  .and(jrcRiverCorridor.unmask(0).not())
  .selfMask()
  .rename("OPTICAL_OUTSIDE_JRC");

var hydroWithoutOptical = hydroRiverCorridor
  .and(opticalWater.unmask(0).not())
  .selfMask()
  .rename("HYDRO_WITHOUT_OPTICAL");

var jrcWithoutOptical = jrcRiverCorridor
  .and(opticalWater.unmask(0).not())
  .selfMask()
  .rename("JRC_WITHOUT_OPTICAL");

print("Optical water pixel count", opticalWater.reduceRegion({
  reducer: ee.Reducer.count(),
  geometry: aoi,
  scale: scaleMeters,
  bestEffort: true,
  maxPixels: 1e9
}).values().get(0));

print("River analysis pixel count", riverAnalysisMask.reduceRegion({
  reducer: ee.Reducer.count(),
  geometry: aoi,
  scale: scaleMeters,
  bestEffort: true,
  maxPixels: 1e9
}).values().get(0));

print("Hydro-constrained river pixel count", hydroConstrainedRiverMask.reduceRegion({
  reducer: ee.Reducer.count(),
  geometry: aoi,
  scale: scaleMeters,
  bestEffort: true,
  maxPixels: 1e9
}).values().get(0));

print("JRC-constrained river pixel count", jrcConstrainedRiverMask.reduceRegion({
  reducer: ee.Reducer.count(),
  geometry: aoi,
  scale: scaleMeters,
  bestEffort: true,
  maxPixels: 1e9
}).values().get(0));

print("JRC occurrence range", jrcOccurrence.reduceRegion({
  reducer: ee.Reducer.percentile([5, 50, 95]),
  geometry: aoi,
  scale: 30,
  bestEffort: true,
  maxPixels: 1e9
}));

print("NDWI on optical water", recent.select("NDWI").updateMask(opticalWater).reduceRegion({
  reducer: ee.Reducer.percentile([5, 50, 95]),
  geometry: aoi,
  scale: scaleMeters,
  bestEffort: true,
  maxPixels: 1e9
}));

print("MNDWI on optical water", recent.select("MNDWI").updateMask(opticalWater).reduceRegion({
  reducer: ee.Reducer.percentile([5, 50, 95]),
  geometry: aoi,
  scale: scaleMeters,
  bestEffort: true,
  maxPixels: 1e9
}));

// ================= WATER-ONLY PROXIES =================
var ndtiWater = recent.select("NDTI")
  .updateMask(riverAnalysisMask)
  .rename("NDTI_WATER");

var tssWater = recent.select("TSS_PROXY")
  .updateMask(riverAnalysisMask)
  .rename("TSS_WATER");

var redGreenWater = recent.select("RED_GREEN")
  .updateMask(riverAnalysisMask)
  .rename("RED_GREEN_WATER");

var waterProxyComposite = ee.Image.cat([
  recent.select("TSS_PROXY").unitScale(20, 220).clamp(0, 1),
  recent.select("NDTI").unitScale(-0.1, 0.25).clamp(0, 1),
  recent.select("RED_GREEN").unitScale(0.7, 2.0).clamp(0, 1)
]).updateMask(riverAnalysisMask).rename(["TSS_RGB", "NDTI_RGB", "RED_GREEN_RGB"]);

// ================= VISUALIZATION =================
var visTrueColor = {bands: ["B4", "B3", "B2"], min: 0.02, max: 0.30};
var visFalseColor = {bands: ["B8", "B4", "B3"], min: 0.03, max: 0.45};
var visObsCount = {min: 0, max: 30, palette: ["2b2b2b", "f7f7f7", "00ff00"]};
var visOccurrence = {min: 0, max: 100, palette: ["0B1F3A", "4FC3F7", "E1F5FE"]};
var visTurbidity = {min: 20, max: 300, palette: ["08306b", "41b6c4", "ffffbf", "fdae61", "d73027"]};
var visNDTI = {min: -0.3, max: 0.4, palette: ["2166AC", "F7F7F7", "B2182B"]};
var visRatio = {min: 0.7, max: 2.0, palette: ["F7FBFF", "FDAE61", "D73027"]};
var visProxyComposite = {bands: ["TSS_RGB", "NDTI_RGB", "RED_GREEN_RGB"], min: 0, max: 1};

// ================= MAP LAYERS =================
Map.addLayer(recent, visTrueColor, "Recent composite (true color)", true);
Map.addLayer(recent, visFalseColor, "Recent composite (false color NIR)", false);
Map.addLayer(obsCount, visObsCount, "QA observation count", false);

Map.addLayer(jrcOccurrence, visOccurrence, "JRC water occurrence", false);
Map.addLayer(jrcRiverSeed, {palette: ["42A5F5"]}, "JRC river seed", false);
Map.addLayer(jrcRiverCorridor, {palette: ["90CAF9"]}, "JRC river corridor", false);
Map.addLayer(hydroRiverSeed, {palette: ["1E88E5"]}, "HydroSHEDS river seed", false);
Map.addLayer(hydroRiverCorridor, {palette: ["64B5F6"]}, "HydroSHEDS river corridor", false);
Map.addLayer(opticalWater, {palette: ["00BFFF"]}, "Optical water mask (primary)", true);
Map.addLayer(riverAnalysisMask, {palette: ["FFD54F"]}, "River analysis mask (optical-first)", true);
Map.addLayer(hydroConstrainedRiverMask, {palette: ["FF7043"]}, "Hydro-constrained river mask", false);
Map.addLayer(jrcConstrainedRiverMask, {palette: ["AB47BC"]}, "JRC-constrained river mask", false);
Map.addLayer(opticalOutsideHydro, {palette: ["FF1744"]}, "QA optical water outside HydroSHEDS corridor", false);
Map.addLayer(opticalOutsideJrc, {palette: ["FF00FF"]}, "QA optical water outside JRC corridor", false);
Map.addLayer(hydroWithoutOptical, {palette: ["FFB300"]}, "QA HydroSHEDS corridor without optical water", false);
Map.addLayer(jrcWithoutOptical, {palette: ["FF7043"]}, "QA JRC corridor without optical water", false);

Map.addLayer(waterProxyComposite, visProxyComposite, "Water proxy composite (TSS / NDTI / red-green)", true);
Map.addLayer(tssWater, visTurbidity, "Water TSS proxy", false);
Map.addLayer(ndtiWater, visNDTI, "Water NDTI", false);
Map.addLayer(redGreenWater, visRatio, "Water red/green ratio", false);

// ================= PANEL =================
var panel = ui.Panel({
  style: {position: "top-right", width: "350px", padding: "8px"}
});

panel.add(ui.Label({
  value: "San Sebastian River Characterization",
  style: {fontWeight: "bold", fontSize: "13px"}
}));
panel.add(ui.Label("Debug order:"));
panel.add(ui.Label("[x] Optical water mask (primary river candidate)"));
panel.add(ui.Label("[x] HydroSHEDS river corridor"));
panel.add(ui.Label("[x] JRC river corridor"));
panel.add(ui.Label("[x] Final river analysis mask (optical-first)"));
panel.add(ui.Label("[x] Water proxy composite"));
panel.add(ui.Label({
  value: "JRC and HydroSHEDS are QA/context layers here, not hard gates. If the optical mask looks right, trust that before adding constraints.",
  style: {fontSize: "11px", color: "B22222"}
}));

Map.add(panel);

print("River characterization script ready. Start with the optical water mask, then compare HydroSHEDS and JRC as QA/context layers.");
