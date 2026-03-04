/*
  FILE: Phase I/gee_scripts/turbidity/san_sebastian_river_mask_visual_simple.js
  PURPOSE: Minimal visual QA script for river-mask inspection at San Sebastian.

  USE THIS FOR
  - quickly checking whether the river mask matches visible channel pixels
  - tuning only a couple of thresholds before any upstream/downstream analysis
*/

// ================= USER SETTINGS =================
var siteName = "San Sebastian Mine (MRDS)";
var siteLon = -87.93002;
var siteLat = 13.6509;

var aoiRadiusMeters = 5000;
var recentMonths = 12;
var cloudMax = 60;
var minObsCount = 3;
var useWetSeasonFilter = true;
var wetSeasonMonths = [5, 6, 7, 8, 9, 10];

// Keep tuning minimal: only these two controls usually matter most.
var opticalOccurrenceThreshold = 0.10;  // fraction of images flagged as water
var ndwiThreshold = 0.05;
var mndwiThreshold = 0.02;
var minConnectedPixels = 8;
var hydroAccThreshold = 15;
var hydroAssistOccurrenceThreshold = 0.06;
var hydroCorridorExpandMeters = 80;

// ================= MAP SETUP =================
Map.setOptions("SATELLITE");
var sitePoint = ee.Geometry.Point([siteLon, siteLat]);
var aoi = sitePoint.buffer(aoiRadiusMeters);
Map.centerObject(sitePoint, 13);
Map.addLayer(sitePoint, {color: "FF0000"}, siteName, true);
Map.addLayer(aoi, {color: "FFFFFF"}, "AOI", false);

// ================= DATE WINDOW =================
var endDate = ee.Date(Date.now());
var startDate = endDate.advance(-recentMonths, "month");
print("Recent window", startDate.format("YYYY-MM-dd"), endDate.format("YYYY-MM-dd"));
print("Occurrence threshold", opticalOccurrenceThreshold);
print("NDWI threshold", ndwiThreshold);
print("MNDWI threshold", mndwiThreshold);
print("Minimum observations", minObsCount);
print("Wet-season month filter", useWetSeasonFilter ? wetSeasonMonths : "OFF");
print("Minimum connected pixels", minConnectedPixels);
print("HydroSHEDS accumulation threshold", hydroAccThreshold);
print("Hydro assist occurrence threshold", hydroAssistOccurrenceThreshold);

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

function addWaterBands(img) {
  var ndwi = img.normalizedDifference(["B3", "B8"]).rename("NDWI");
  var mndwi = img.normalizedDifference(["B3", "B11"]).rename("MNDWI");
  return img.addBands([ndwi, mndwi]);
}

function addMonth(img) {
  return img.set("month", ee.Date(img.get("system:time_start")).get("month"));
}

var s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  .filterBounds(aoi)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", cloudMax))
  .map(maskS2)
  .map(scaleS2)
  .map(addMonth)
  .map(addWaterBands);

if (useWetSeasonFilter) {
  s2 = s2.filter(ee.Filter.inList("month", wetSeasonMonths));
}

print("Recent image count", s2.size());

var recent = s2.median().clip(aoi);
var obsCount = s2.select("B4").count().clip(aoi).rename("OBS_COUNT");
var enoughObs = obsCount.gte(minObsCount);

var opticalInstant = recent.select("NDWI").gt(ndwiThreshold)
  .or(recent.select("MNDWI").gt(mndwiThreshold))
  .updateMask(enoughObs)
  .selfMask()
  .rename("OPTICAL_INSTANT");

var opticalFraction = s2.map(function(img) {
  var w = img.select("NDWI").gt(ndwiThreshold).or(img.select("MNDWI").gt(mndwiThreshold));
  return w.rename("WATER_BIN").toFloat();
}).mean().clip(aoi).rename("OPTICAL_WATER_FRACTION");

var opticalOccurrence = opticalFraction.gte(opticalOccurrenceThreshold)
  .updateMask(enoughObs)
  .selfMask()
  .rename("OPTICAL_OCCURRENCE");

var rawRiverMask = opticalOccurrence.or(opticalInstant).selfMask().rename("RIVER_MASK_RAW");

var connectedMask = rawRiverMask.connectedPixelCount(100, true)
  .gte(minConnectedPixels)
  .selfMask()
  .rename("CONNECTED_MASK");

var filteredRiverMask = rawRiverMask
  .updateMask(connectedMask)
  .focal_max(20, "circle", "meters")
  .focal_min(20, "circle", "meters")
  .selfMask()
  .rename("RIVER_MASK_FILTERED");

var hydroAcc = ee.Image("WWF/HydroSHEDS/15ACC").clip(aoi).rename("HYDRO_ACC");
var hydroCorridor = hydroAcc.gte(hydroAccThreshold)
  .focal_max(hydroCorridorExpandMeters, "circle", "meters")
  .selfMask()
  .rename("HYDRO_CORRIDOR");

var hydroAssistMask = hydroCorridor
  .and(opticalFraction.gte(hydroAssistOccurrenceThreshold))
  .selfMask()
  .rename("HYDRO_ASSIST_MASK");

var riverMask = filteredRiverMask.or(hydroAssistMask).selfMask().rename("RIVER_MASK");

print("Instant optical pixel count", opticalInstant.reduceRegion({
  reducer: ee.Reducer.count(),
  geometry: aoi,
  scale: 20,
  bestEffort: true,
  maxPixels: 1e9
}).values().get(0));

print("Occurrence optical pixel count", opticalOccurrence.reduceRegion({
  reducer: ee.Reducer.count(),
  geometry: aoi,
  scale: 20,
  bestEffort: true,
  maxPixels: 1e9
}).values().get(0));

print("Filtered river-mask pixel count", filteredRiverMask.reduceRegion({
  reducer: ee.Reducer.count(),
  geometry: aoi,
  scale: 20,
  bestEffort: true,
  maxPixels: 1e9
}).values().get(0));

print("Final river-mask pixel count", riverMask.reduceRegion({
  reducer: ee.Reducer.count(),
  geometry: aoi,
  scale: 20,
  bestEffort: true,
  maxPixels: 1e9
}).values().get(0));

// ================= VISUALIZATION =================
var visTrueColor = {bands: ["B4", "B3", "B2"], min: 0.02, max: 0.30};
var visFalseColor = {bands: ["B8", "B4", "B3"], min: 0.03, max: 0.45};
var visObs = {min: 0, max: 30, palette: ["2b2b2b", "f7f7f7", "00ff00"]};
var visFraction = {min: 0, max: 0.35, palette: ["1B1B1B", "00ACC1", "FFF176"]};

Map.addLayer(recent, visTrueColor, "Recent composite (true color)", true);
Map.addLayer(recent, visFalseColor, "Recent composite (false color)", false);
Map.addLayer(obsCount, visObs, "QA observation count", false);
Map.addLayer(opticalFraction, visFraction, "Optical water occurrence fraction", true);
Map.addLayer(opticalInstant, {palette: ["00B8D4"]}, "Optical instant mask", false);
Map.addLayer(opticalOccurrence, {palette: ["26C6DA"]}, "Optical occurrence mask", false);
Map.addLayer(rawRiverMask, {palette: ["FFF176"]}, "River mask raw", false);
Map.addLayer(filteredRiverMask, {palette: ["FFB300"]}, "River mask filtered", false);
Map.addLayer(hydroCorridor, {palette: ["64B5F6"]}, "HydroSHEDS corridor", false);
Map.addLayer(hydroAssistMask, {palette: ["1E88E5"]}, "Hydro-assisted mask", false);
Map.addLayer(riverMask, {palette: ["FFD54F"]}, "River mask (visual QA)", true);

var panel = ui.Panel({style: {position: "top-right", width: "340px", padding: "8px"}});
panel.add(ui.Label({
  value: "River Mask Visual QA (Simple)",
  style: {fontWeight: "bold", fontSize: "13px"}
}));
panel.add(ui.Label("Tune only these first:"));
panel.add(ui.Label("1) opticalOccurrenceThreshold"));
panel.add(ui.Label("2) ndwiThreshold / mndwiThreshold"));
panel.add(ui.Label("3) minConnectedPixels"));
panel.add(ui.Label({
  value: "Goal: yellow mask should cover visible channel portions without flooding adjacent land.",
  style: {fontSize: "11px", color: "B22222"}
}));
Map.add(panel);
