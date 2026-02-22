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

var aoiRadiusMeters = 3000;   // narrowed scope for faster, site-focused analysis
var innerFocusMeters = 1000;
var recentMonths = 12;      // "current" window; keeps at least one dry-season pass most of the year
var baselineYears = 5;      // prior multi-year baseline for anomaly comparison
var cloudMax = 60;
var seasonMonths = [11, 12, 1, 2, 3, 4]; // drier season often better for turbidity visibility
var useSeasonFilterForCurrentMaps = true;
var enableTimeSeriesCharts = false; // set true when you want the heavier monthly analysis

// Anomaly / significance filtering (baseline-vs-current)
var sigmaMin = 0.02;           // floor for std dev to avoid unstable z-scores
var zSigThreshold = 2.0;       // "significant" anomaly cutoff
var anomalyHeatThreshold = 0.55; // combined anomaly heat threshold (0-1)

// Upstream/downstream proxy corridor settings (user-tunable)
// Azimuth is downstream direction in degrees clockwise from north.
// Set from local river direction after visual inspection.
var downstreamAzimuthDeg = 110;
var riverCorridorExpandMeters = 120;

// Time-series settings
var seriesYears = 3;
var useSeasonFilterForTimeSeries = false; // monthly tracking usually should include all months

// ================= MAP SETUP =================
Map.setOptions("SATELLITE");

var minePoint = ee.Geometry.Point([mineLon, mineLat]);
var aoi = minePoint.buffer(aoiRadiusMeters);
var innerFocus = minePoint.buffer(innerFocusMeters);

Map.centerObject(minePoint, 12);
Map.addLayer(aoi, {color: "FFFFFF"}, "AOI (3 km)", true);
Map.addLayer(innerFocus, {color: "FFAA00"}, "Inner focus (1 km)", false);
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

function buildS2Collection(start, end, applySeasonFilter) {
  var col = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
    .filterBounds(aoi)
    .filterDate(start, end)
    .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", cloudMax))
    .map(maskS2)
    .map(scaleS2)
    .map(addMonth);

  if (applySeasonFilter) {
    col = col.filter(ee.Filter.inList("month", seasonMonths));
  }

  return col.map(addIndices);
}

function compositeWithCounts(col, label) {
  print(label + " image count", col.size());
  return col.median().clip(aoi);
}

var recentCol = buildS2Collection(recentStart, recentEnd, useSeasonFilterForCurrentMaps);
var baselineCol = buildS2Collection(baselineStart, baselineEnd, true);

var recent = compositeWithCounts(recentCol, "Recent");
var baseline = compositeWithCounts(baselineCol, "Baseline");
var recentObsCount = recentCol.select("B4").count().clip(aoi).rename("RECENT_OBS_COUNT");
var baselineObsCount = baselineCol.select("B4").count().clip(aoi).rename("BASELINE_OBS_COUNT");

var baselineMean = baselineCol.select([
  "NDVI", "SAVI", "NDMI", "BSI", "NDTI", "IOI", "FERROUS", "CLAY", "TSS_PROXY"
]).mean().clip(aoi);

var baselineStd = baselineCol.select([
  "NDVI", "SAVI", "NDMI", "BSI", "NDTI", "IOI", "FERROUS", "CLAY", "TSS_PROXY"
]).reduce(ee.Reducer.stdDev()).clip(aoi);

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

// ================= HYDROLOGY-ORIENTED PROXY CORRIDORS =================
// Proxy approach: split a buffered water corridor into upstream/downstream halves
// using a user-defined downstream azimuth. This is not full flow-routing.
var riverCorridor = permanentWater.focal_max(riverCorridorExpandMeters, "circle", "meters")
  .clip(aoi)
  .selfMask()
  .rename("RIVER_CORRIDOR");

var azimuthRad = downstreamAzimuthDeg * Math.PI / 180;
var ux = Math.sin(azimuthRad); // east-west component
var uy = Math.cos(azimuthRad); // north-south component

var lonLat = ee.Image.pixelLonLat();
var signedAlong = lonLat.select("longitude").subtract(mineLon).multiply(ux)
  .add(lonLat.select("latitude").subtract(mineLat).multiply(uy))
  .rename("SIGNED_ALONG");

var downstreamHalf = signedAlong.gt(0);
var upstreamHalf = signedAlong.lte(0);

var downstreamCorridor = riverCorridor.and(downstreamHalf).selfMask().rename("DOWNSTREAM_CORRIDOR");
var upstreamCorridor = riverCorridor.and(upstreamHalf).selfMask().rename("UPSTREAM_CORRIDOR");

var downstreamWaterNow = stableWaterMask.and(downstreamCorridor).selfMask();
var upstreamWaterNow = stableWaterMask.and(upstreamCorridor).selfMask();

// ================= SCREENING SCORES (RED = HIGHER CONCERN) =================
function norm(img, min, max) {
  return img.unitScale(min, max).clamp(0, 1).rename(img.bandNames());
}

function invertNorm(img, min, max) {
  return ee.Image(1).subtract(norm(img, min, max)).rename(img.bandNames());
}

function zPos(diffImg, stdImg, outName) {
  var safeStd = stdImg.max(sigmaMin);
  return diffImg.divide(safeStd).max(0).rename(outName);
}

var landMask = stableWaterMask.not();
var waterMaskImg = stableWaterMask;

// ================= BASELINE ANOMALIES (SIGNIFICANCE-FOCUSED) =================
// Positive values are "more concern" relative to baseline.
var bsiAnom = recent.select("BSI").subtract(baselineMean.select("BSI")).rename("BSI_ANOM");
var ndviDrop = baselineMean.select("NDVI").subtract(recent.select("NDVI")).rename("NDVI_DROP");
var ndmiDrop = baselineMean.select("NDMI").subtract(recent.select("NDMI")).rename("NDMI_DROP");
var ioiAnom = recent.select("IOI").subtract(baselineMean.select("IOI")).rename("IOI_ANOM");
var ferrousAnom = recent.select("FERROUS").subtract(baselineMean.select("FERROUS")).rename("FERROUS_ANOM");
var clayAnom = recent.select("CLAY").subtract(baselineMean.select("CLAY")).rename("CLAY_ANOM");
var ndtiAnom = recent.select("NDTI").subtract(baselineMean.select("NDTI")).rename("NDTI_ANOM");
var tssAnomRaw = recent.select("TSS_PROXY").subtract(baselineMean.select("TSS_PROXY")).rename("TSS_ANOM_BASE");

var bsiZ = zPos(bsiAnom, baselineStd.select("BSI_stdDev"), "BSI_ZP").updateMask(landMask);
var ndviDropZ = zPos(ndviDrop, baselineStd.select("NDVI_stdDev"), "NDVI_DROP_ZP").updateMask(landMask);
var ndmiDropZ = zPos(ndmiDrop, baselineStd.select("NDMI_stdDev"), "NDMI_DROP_ZP").updateMask(landMask);
var ioiZ = zPos(ioiAnom, baselineStd.select("IOI_stdDev"), "IOI_ZP").updateMask(landMask);
var ferrousZ = zPos(ferrousAnom, baselineStd.select("FERROUS_stdDev"), "FERROUS_ZP").updateMask(landMask);
var clayZ = zPos(clayAnom, baselineStd.select("CLAY_stdDev"), "CLAY_ZP").updateMask(landMask);
var ndtiWaterZ = zPos(ndtiAnom, baselineStd.select("NDTI_stdDev"), "NDTI_W_ZP").updateMask(waterMaskImg);
var tssWaterZ = zPos(tssAnomRaw, baselineStd.select("TSS_PROXY_stdDev"), "TSS_W_ZP").updateMask(waterMaskImg);

var ndviDropMask = ndviDropZ.gte(zSigThreshold).selfMask().rename("NDVI_DROP_SIG");
var bsiRiseMask = bsiZ.gte(zSigThreshold).selfMask().rename("BSI_RISE_SIG");
var ioiRiseMask = ioiZ.gte(zSigThreshold).selfMask().rename("IOI_RISE_SIG");
var turbidityRiseMask = tssWaterZ.gte(zSigThreshold).selfMask().rename("TSS_RISE_SIG");
var waterNDTIMask = ndtiWaterZ.gte(zSigThreshold).selfMask().rename("NDTI_W_SIG");

var landAnomalyHeat = ee.ImageCollection.fromImages([
  norm(bsiZ, 0, 4),
  norm(ndviDropZ, 0, 4),
  norm(ndmiDropZ, 0, 4),
  norm(ioiZ, 0, 4),
  norm(ferrousZ, 0, 4),
  norm(clayZ, 0, 4)
]).mean().updateMask(landMask).rename("LAND_ANOM_HEAT");

var waterAnomalyHeat = ee.ImageCollection.fromImages([
  norm(tssWaterZ, 0, 4),
  norm(ndtiWaterZ, 0, 4)
]).mean().updateMask(waterMaskImg).rename("WATER_ANOM_HEAT");

var anomalyHeat = landAnomalyHeat.unmask(0).max(waterAnomalyHeat.unmask(0))
  .rename("ANOMALY_HEAT");
var anomalySigMask = anomalyHeat.gte(anomalyHeatThreshold).selfMask().rename("ANOMALY_SIG_MASK");

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
var visObsCount = {min: 0, max: 30, palette: ["2b2b2b", "f7f7f7", "00ff00"]};
var visZ = {min: 0, max: 4, palette: ["f7f7f7", "fdae61", "d73027", "7f0000"]};
var visAnomHeat = {min: 0, max: 1, palette: ["fff7ec", "fdd49e", "fc8d59", "d7301f", "7f0000"]};

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

// Baseline-vs-current anomaly layers (focus on significant deviations)
Map.addLayer(bsiAnom.updateMask(landMask), {min: -0.2, max: 0.3, palette: ["2166ac", "f7f7f7", "b2182b"]}, "BSI anomaly vs baseline", false);
Map.addLayer(ndviDrop.updateMask(landMask), {min: -0.3, max: 0.4, palette: ["2166ac", "f7f7f7", "b2182b"]}, "NDVI drop vs baseline (red=drop)", false);
Map.addLayer(bsiZ, visZ, "BSI anomaly significance (z+)", false);
Map.addLayer(ndviDropZ, visZ, "NDVI-drop significance (z+)", false);
Map.addLayer(ioiZ, visZ, "IOI anomaly significance (z+)", false);
Map.addLayer(ndtiWaterZ, visZ, "Water NDTI anomaly significance (z+)", false);
Map.addLayer(tssWaterZ, visZ, "Turbidity anomaly significance (z+)", true);
Map.addLayer(landAnomalyHeat, visAnomHeat, "Land anomaly heatmap", true);
Map.addLayer(waterAnomalyHeat, visAnomHeat, "Water anomaly heatmap", true);
Map.addLayer(anomalyHeat, visAnomHeat, "Combined anomaly heatmap", true);
Map.addLayer(anomalySigMask, {palette: ["FF0000"]}, "Significant anomaly mask (filtered)", true);
Map.addLayer(ndviDropMask, {palette: ["FF6F00"]}, "Mask: significant NDVI drop", false);
Map.addLayer(bsiRiseMask, {palette: ["FF0000"]}, "Mask: significant BSI rise", false);
Map.addLayer(ioiRiseMask, {palette: ["C2185B"]}, "Mask: significant IOI rise", false);

// Water and turbidity
Map.addLayer(stableWaterMask.selfMask(), {palette: ["00BFFF"]}, "Stable water mask", false);
Map.addLayer(riverCorridor, {palette: ["80DEEA"]}, "River corridor (proxy buffer)", false);
Map.addLayer(upstreamCorridor, {palette: ["1E88E5"]}, "Upstream proxy corridor", false);
Map.addLayer(downstreamCorridor, {palette: ["E53935"]}, "Downstream proxy corridor", false);
Map.addLayer(recentTurbidity, visTurbidity, "Turbidity/TSS proxy (water)", true);
Map.addLayer(turbidityAnomaly, visTurbidityAnom, "Turbidity anomaly vs baseline (red=increased)", true);
Map.addLayer(recentNDTIWater, {min: -0.1, max: 0.3, palette: ["2166ac", "f7f7f7", "d73027"]}, "Water NDTI (discoloration)", false);
Map.addLayer(turbidityRiseMask, {palette: ["FF1744"]}, "Mask: significant turbidity rise", true);
Map.addLayer(waterNDTIMask, {palette: ["D50000"]}, "Mask: significant water NDTI rise", false);
Map.addLayer(upstreamWaterNow, {palette: ["2962FF"]}, "Upstream water (current)", false);
Map.addLayer(downstreamWaterNow, {palette: ["FF1744"]}, "Downstream water (current)", false);

// Screening layers (red-forward for quick issue scan)
Map.addLayer(landStress, visStress, "Land mine-impact screening (red=high)", true);
Map.addLayer(waterStress, visStress, "Water contamination screening (red=high)", true);
Map.addLayer(combinedScreen, visStress, "Combined environmental screening (red=high)", true);
Map.addLayer(hotSpots, {palette: ["FF0000"]}, "Hotspots > 0.65", true);
Map.addLayer(recentObsCount, visObsCount, "QA: recent observation count", false);
Map.addLayer(baselineObsCount, visObsCount, "QA: baseline observation count", false);

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
panel.add(ui.Label("[x] Baseline-vs-current anomaly heatmap (filtered)"));
panel.add(ui.Label("[x] Significant anomaly masks (z-threshold)"));
panel.add(ui.Label("[x] Upstream/downstream proxy river corridor split"));
panel.add(ui.Label("[x] Monthly upstream vs downstream time-series charts"));
panel.add(ui.Label("[x] QA observation-count layers"));

panel.add(ui.Label({
  value: "Red means higher screening concern (disturbance/stress/turbidity proxy), not confirmed sulfide/nickel concentration.",
  style: {fontSize: "11px", color: "B22222"}
}));
panel.add(ui.Label({
  value: "Anomaly filter uses baseline-vs-current significance (z >= " + zSigThreshold + "), combined heat threshold = " + anomalyHeatThreshold + ".",
  style: {fontSize: "11px"}
}));

Map.add(panel);

print("NOTE: Sulfide, nickel, gold, and dissolved metals require water/soil chemistry sampling.");
print("This script is for remote sensing triage and hotspot targeting.");
print("Anomaly thresholds", ee.Dictionary({
  zSigThreshold: zSigThreshold,
  anomalyHeatThreshold: anomalyHeatThreshold,
  sigmaMin: sigmaMin
}));

// ================= MONTHLY UPSTREAM/DOWNSTREAM TIME SERIES =================
function maskedMean(img, bandName, maskImg) {
  var stats = img.select(bandName).updateMask(maskImg).reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi,
    scale: 10,
    bestEffort: true,
    maxPixels: 1e9
  });
  return stats.get(bandName);
}

var monthlyStats = ee.FeatureCollection([]);
if (enableTimeSeriesCharts) {
  var seriesStart = ee.Date(now.advance(-seriesYears, "year").format("YYYY-MM-01"));
  var monthCount = ee.Number(now.difference(seriesStart, "month")).floor();
  var monthOffsets = ee.List.sequence(0, monthCount.subtract(1));

  monthlyStats = ee.FeatureCollection(monthOffsets.map(function(m) {
    m = ee.Number(m);
    var start = seriesStart.advance(m, "month");
    var end = start.advance(1, "month");
    var col = buildS2Collection(start, end, useSeasonFilterForTimeSeries);
    var count = col.size();

    var precip = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY")
      .filterBounds(aoi)
      .filterDate(start, end)
      .sum();
    var rainMm = precip.reduceRegion({
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
      down_tss: null,
      up_ndti_w: null,
      down_ndti_w: null,
      up_minus_down_tss: null,
      down_minus_up_tss: null,
      up_obs_px: 0,
      down_obs_px: 0
    });

    var populatedFeature = ee.Feature(ee.Algorithms.If(count.gt(0), (function() {
      var img = ee.Image(col.median()).clip(aoi);
      var w = waterMask(img);
      var upMask = w.and(upstreamCorridor);
      var downMask = w.and(downstreamCorridor);

      var upTss = maskedMean(img, "TSS_PROXY", upMask);
      var downTss = maskedMean(img, "TSS_PROXY", downMask);
      var upNDTI = maskedMean(img, "NDTI", upMask);
      var downNDTI = maskedMean(img, "NDTI", downMask);

      var upObsPx = ee.Number(upMask.reduceRegion({
        reducer: ee.Reducer.count(),
        geometry: aoi,
        scale: 10,
        bestEffort: true,
        maxPixels: 1e9
      }).get("NDWI", 0));

      var downObsPx = ee.Number(downMask.reduceRegion({
        reducer: ee.Reducer.count(),
        geometry: aoi,
        scale: 10,
        bestEffort: true,
        maxPixels: 1e9
      }).get("NDWI", 0));

      var upMinusDown = ee.Algorithms.If(
        ee.Algorithms.IsEqual(upTss, null),
        null,
        ee.Algorithms.If(ee.Algorithms.IsEqual(downTss, null), null,
          ee.Number(upTss).subtract(ee.Number(downTss)))
      );

      var downMinusUp = ee.Algorithms.If(
        ee.Algorithms.IsEqual(upTss, null),
        null,
        ee.Algorithms.If(ee.Algorithms.IsEqual(downTss, null), null,
          ee.Number(downTss).subtract(ee.Number(upTss)))
      );

      return ee.Feature(null, {
        "system:time_start": start.millis(),
        date: start.format("YYYY-MM"),
        img_count: count,
        rain_mm: rainMm,
        up_tss: upTss,
        down_tss: downTss,
        up_ndti_w: upNDTI,
        down_ndti_w: downNDTI,
        up_minus_down_tss: upMinusDown,
        down_minus_up_tss: downMinusUp,
        up_obs_px: upObsPx,
        down_obs_px: downObsPx
      });
    })(), emptyFeature));

    return populatedFeature;
  })).sort("system:time_start");
}

if (enableTimeSeriesCharts) {
  print("Monthly upstream/downstream stats table", monthlyStats.limit(12));

  var tssChart = ui.Chart.feature.byFeature(monthlyStats, "date", ["up_tss", "down_tss"])
    .setChartType("LineChart")
    .setOptions({
      title: "Monthly Water TSS Proxy (Upstream vs Downstream proxy corridors)",
      hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
      vAxis: {title: "TSS proxy (relative)"},
      lineWidth: 2,
      pointSize: 3,
      series: {
        0: {color: "#2962FF"},
        1: {color: "#D50000"}
      }
    });
  print(tssChart);

  var tssDiffChart = ui.Chart.feature.byFeature(monthlyStats, "date", ["down_minus_up_tss"])
    .setChartType("ColumnChart")
    .setOptions({
      title: "Monthly TSS Proxy Difference (Downstream - Upstream)",
      hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
      vAxis: {title: "Difference (relative)"},
      colors: ["#B71C1C"]
    });
  print(tssDiffChart);

  var ndtiChart = ui.Chart.feature.byFeature(monthlyStats, "date", ["up_ndti_w", "down_ndti_w"])
    .setChartType("LineChart")
    .setOptions({
      title: "Monthly Water NDTI (Discoloration Proxy)",
      hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
      vAxis: {title: "NDTI"},
      lineWidth: 2,
      pointSize: 3,
      series: {
        0: {color: "#1565C0"},
        1: {color: "#C62828"}
      }
    });
  print(ndtiChart);

  var rainChart = ui.Chart.feature.byFeature(monthlyStats, "date", ["rain_mm"])
    .setChartType("ColumnChart")
    .setOptions({
      title: "Monthly CHIRPS Rainfall (AOI mean, mm)",
      hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
      vAxis: {title: "Rainfall (mm)"},
      colors: ["#26A69A"]
    });
  print(rainChart);

  var qaChart = ui.Chart.feature.byFeature(monthlyStats, "date", ["img_count", "up_obs_px", "down_obs_px"])
    .setChartType("LineChart")
    .setOptions({
      title: "QA: Monthly image count and water pixels used",
      hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
      vAxis: {title: "Count"},
      lineWidth: 2,
      pointSize: 2,
      series: {
        0: {color: "#424242"},
        1: {color: "#1E88E5"},
        2: {color: "#E53935"}
      }
    });
  print(qaChart);
} else {
  print("Monthly charts disabled (enableTimeSeriesCharts=false) for lighter runs.");
}

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
