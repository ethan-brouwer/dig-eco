/*
  FILE: Phase I/gee_scripts/turbidity/san_sebastian_ssrivers_buffer_indices.js
  PURPOSE: Run river-only index screening inside a manual SSrivers corridor.

  REQUIREMENT
  - In GEE Imports, provide `SSrivers` as a Geometry, Feature, or FeatureCollection.
*/

// ================= USER SETTINGS =================
var siteName = "San Sebastian Mine (MRDS)";
var siteLon = -87.93002;
var siteLat = 13.6509;

var recentMonths = 12;
var cloudMax = 60;
var useWetSeasonFilter = true;
var wetSeasonMonths = [5, 6, 7, 8, 9, 10];

var corridorBufferMeters = 120;
var minObsCount = 3;
var ndwiThreshold = 0.05;
var mndwiThreshold = 0.02;
var opticalOccurrenceThreshold = 0.10;

// ================= MAP SETUP =================
Map.setOptions("SATELLITE");
var sitePoint = ee.Geometry.Point([siteLon, siteLat]);
Map.centerObject(sitePoint, 12);
Map.addLayer(sitePoint, {color: "FF0000"}, siteName, true);

// ================= IMPORTED GEOMETRY (SSrivers) =================
var ssType = ee.String(ee.Algorithms.ObjectType(SSrivers));
var ssGeom = ee.Geometry(ee.Algorithms.If(
  ssType.compareTo("FeatureCollection").eq(0),
  ee.FeatureCollection(SSrivers).geometry(),
  ee.Algorithms.If(
    ssType.compareTo("Feature").eq(0),
    ee.Feature(SSrivers).geometry(),
    ee.Geometry(SSrivers)
  )
));

var riverCorridor = ssGeom.buffer(corridorBufferMeters);
var aoi = riverCorridor.bounds().buffer(1000);

Map.addLayer(ssGeom, {color: "00E5FF"}, "SSrivers line/geometry", true);
Map.addLayer(riverCorridor, {color: "FFD54F"}, "Buffered river corridor", true);
Map.addLayer(aoi, {color: "FFFFFF"}, "AOI", false);

// ================= DATE WINDOW =================
var endDate = ee.Date(Date.now());
var startDate = endDate.advance(-recentMonths, "month");
print("Recent window", startDate.format("YYYY-MM-dd"), endDate.format("YYYY-MM-dd"));
print("Wet-season month filter", useWetSeasonFilter ? wetSeasonMonths : "OFF");
print("Corridor buffer (m)", corridorBufferMeters);
print("NDWI threshold", ndwiThreshold);
print("MNDWI threshold", mndwiThreshold);
print("Occurrence threshold", opticalOccurrenceThreshold);

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

function addMonth(img) {
  return img.set("month", ee.Date(img.get("system:time_start")).get("month"));
}

function addIndices(img) {
  var green = img.select("B3");
  var red = img.select("B4");
  var ndwi = img.normalizedDifference(["B3", "B8"]).rename("NDWI");
  var mndwi = img.normalizedDifference(["B3", "B11"]).rename("MNDWI");
  var ndti = img.normalizedDifference(["B4", "B3"]).rename("NDTI");
  var tssProxy = red.multiply(1000).rename("TSS_PROXY");
  var redGreen = red.divide(green).rename("RED_GREEN");
  return img.addBands([ndwi, mndwi, ndti, tssProxy, redGreen]);
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

print("Recent image count", s2.size());

var recent = s2.median().clip(aoi);
var obsCount = s2.select("B4").count().clip(aoi).rename("OBS_COUNT");
var enoughObs = obsCount.gte(minObsCount);

// Optical river detection only inside your manual corridor.
var corridorMask = ee.Image.constant(1).clip(riverCorridor).selfMask().rename("CORRIDOR_MASK");
var opticalInstant = recent.select("NDWI").gt(ndwiThreshold)
  .or(recent.select("MNDWI").gt(mndwiThreshold))
  .updateMask(enoughObs)
  .updateMask(corridorMask)
  .selfMask()
  .rename("OPTICAL_INSTANT");

var opticalFraction = s2.map(function(img) {
  var w = img.select("NDWI").gt(ndwiThreshold).or(img.select("MNDWI").gt(mndwiThreshold));
  return w.rename("WATER_BIN").toFloat();
}).mean().clip(aoi).rename("OPTICAL_WATER_FRACTION");

var opticalOccurrence = opticalFraction.gte(opticalOccurrenceThreshold)
  .updateMask(enoughObs)
  .updateMask(corridorMask)
  .selfMask()
  .rename("OPTICAL_OCCURRENCE");

var riverMask = opticalOccurrence.or(opticalInstant).selfMask().rename("RIVER_MASK");

// ================= WATER-ONLY INDEX LAYERS =================
var ndtiWater = recent.select("NDTI").updateMask(riverMask).rename("NDTI_WATER");
var tssWater = recent.select("TSS_PROXY").updateMask(riverMask).rename("TSS_WATER");
var redGreenWater = recent.select("RED_GREEN").updateMask(riverMask).rename("RED_GREEN_WATER");

var waterProxyComposite = ee.Image.cat([
  recent.select("TSS_PROXY").unitScale(20, 220).clamp(0, 1),
  recent.select("NDTI").unitScale(-0.1, 0.25).clamp(0, 1),
  recent.select("RED_GREEN").unitScale(0.7, 2.0).clamp(0, 1)
]).updateMask(riverMask).rename(["TSS_RGB", "NDTI_RGB", "RED_GREEN_RGB"]);

print("River-mask pixel count", riverMask.reduceRegion({
  reducer: ee.Reducer.count(),
  geometry: riverCorridor,
  scale: 20,
  bestEffort: true,
  maxPixels: 1e9
}).values().get(0));

print("Water index percentiles in corridor", recent.select(["TSS_PROXY", "NDTI", "RED_GREEN"])
  .updateMask(riverMask)
  .reduceRegion({
    reducer: ee.Reducer.percentile([5, 50, 95]),
    geometry: riverCorridor,
    scale: 20,
    bestEffort: true,
    maxPixels: 1e9
  })
);

// ================= MAP LAYERS =================
var visTrueColor = {bands: ["B4", "B3", "B2"], min: 0.02, max: 0.30};
var visFalseColor = {bands: ["B8", "B4", "B3"], min: 0.03, max: 0.45};
var visObs = {min: 0, max: 30, palette: ["2b2b2b", "f7f7f7", "00ff00"]};
var visFraction = {min: 0, max: 0.35, palette: ["1B1B1B", "00ACC1", "FFF176"]};
var visNDTI = {min: -0.3, max: 0.4, palette: ["2166AC", "F7F7F7", "B2182B"]};
var visTSS = {min: 20, max: 300, palette: ["08306b", "41b6c4", "ffffbf", "fdae61", "d73027"]};
var visRatio = {min: 0.7, max: 2.0, palette: ["F7FBFF", "FDAE61", "D73027"]};
var visComposite = {bands: ["TSS_RGB", "NDTI_RGB", "RED_GREEN_RGB"], min: 0, max: 1};

Map.addLayer(recent, visTrueColor, "Recent composite (true color)", true);
Map.addLayer(recent, visFalseColor, "Recent composite (false color)", false);
Map.addLayer(obsCount, visObs, "QA observation count", false);
Map.addLayer(corridorMask, {palette: ["FFF59D"]}, "Corridor mask", false);
Map.addLayer(opticalFraction.updateMask(corridorMask), visFraction, "Optical water occurrence (corridor)", true);
Map.addLayer(opticalInstant, {palette: ["00B8D4"]}, "Optical instant (corridor)", false);
Map.addLayer(opticalOccurrence, {palette: ["26C6DA"]}, "Optical occurrence (corridor)", false);
Map.addLayer(riverMask, {palette: ["FFD54F"]}, "River mask (from SSrivers buffer)", true);
Map.addLayer(waterProxyComposite, visComposite, "Water proxy composite", true);
Map.addLayer(tssWater, visTSS, "Water TSS proxy", false);
Map.addLayer(ndtiWater, visNDTI, "Water NDTI", false);
Map.addLayer(redGreenWater, visRatio, "Water red/green ratio", false);
