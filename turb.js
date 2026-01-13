// Turbidity proxy workflow (sediment/optical only, not chemical contamination)
// Results indicate relative turbidity changes; seasonal hydrology adds uncertainty.

// === AOI ===
var aoi = ee.FeatureCollection("FAO/GAUL/2015/level1")
  .filter(ee.Filter.eq("ADM0_NAME", "El Salvador"))
  .filter(ee.Filter.eq("ADM1_NAME", "Cabanas"))
  .geometry();
Map.centerObject(aoi, 9);

// === Parameters ===
var years = [2018, 2019, 2020, 2021, 2022, 2023];
var dryStart = "-11-01";
var dryEnd = "-04-30";
var cloudMax = 70;
var ndwiThresh = 0.15;
var mndwiThresh = 0.2;
var useSentinel1 = false;

// === Cloud and shadow masking (strict) ===
function maskS2Strict(img) {
  var scl = img.select("SCL");
  var qa = img.select("QA60");
  var cloud = qa.bitwiseAnd(1 << 10).neq(0).or(qa.bitwiseAnd(1 << 11).neq(0));
  var shadow = scl.eq(3);
  var cirrus = scl.eq(10);
  var snow = scl.eq(11);
  var mask = cloud.or(shadow).or(cirrus).or(snow);
  return img.updateMask(mask.not());
}

// === Turbidity features ===
function addTurbidityBands(img) {
  var red = img.select("B4").rename("RED");
  var green = img.select("B3");
  var nir = img.select("B8");
  var redGreen = red.divide(green).rename("RED_GREEN");
  var redNir = red.divide(nir).rename("RED_NIR");
  var tssProxy = red.multiply(1000).rename("TSS_PROXY"); // relative proxy
  var ndwi = img.normalizedDifference(["B3", "B8"]).rename("NDWI");
  var mndwi = img.normalizedDifference(["B3", "B11"]).rename("MNDWI");
  return img.addBands([red, redGreen, redNir, tssProxy, ndwi, mndwi]);
}

function seasonalComposite(year) {
  var start = ee.Date(year + dryStart);
  var end = ee.Date((year + 1) + dryEnd);
  return ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
    .filterBounds(aoi)
    .filterDate(start, end)
    .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", cloudMax))
    .map(maskS2Strict)
    .map(addTurbidityBands)
    .median()
    .clip(aoi);
}

// === Water mask ===
var jrcWater = ee.Image("JRC/GSW1_4/Occurrence");
var permanentWater = jrcWater.gte(70);

function waterMask(img) {
  var ndwi = img.select("NDWI");
  var mndwi = img.select("MNDWI");
  var water = ndwi.gt(ndwiThresh).or(mndwi.gt(mndwiThresh));
  return water.updateMask(permanentWater);
}

// === Optional Sentinel-1 for stable water extent ===
function s1WaterStable(start, end) {
  var s1 = ee.ImageCollection("COPERNICUS/S1_GRD")
    .filterBounds(aoi)
    .filterDate(start, end)
    .filter(ee.Filter.eq("instrumentMode", "IW"))
    .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VV"))
    .select("VV")
    .median()
    .clip(aoi);
  return s1.lt(-17); // simple low-backscatter water proxy
}

// === Build annual turbidity composites ===
var annual = years.map(function(y) {
  var img = seasonalComposite(y);
  var wMask = waterMask(img);
  if (useSentinel1) {
    var s1 = s1WaterStable(
      ee.Date(y + dryStart),
      ee.Date((y + 1) + dryEnd)
    );
    wMask = wMask.and(s1);
  }
  var turb = img.select(["RED", "RED_GREEN", "RED_NIR", "TSS_PROXY"])
    .updateMask(wMask)
    .rename(["RED_" + y, "RED_GREEN_" + y, "RED_NIR_" + y, "TSS_" + y]);
  return turb.set("year", y);
});

var annualCol = ee.ImageCollection.fromImages(annual);

// === Normalize turbidity for inter-annual comparison ===
var tssStack = annualCol.select("TSS_.*");
var tssMin = tssStack.reduce(ee.Reducer.min());
var tssMax = tssStack.reduce(ee.Reducer.max());
var tssNorm = tssStack.subtract(tssMin).divide(tssMax.subtract(tssMin))
  .rename(tssStack.bandNames().map(function(b) {
    return ee.String(b).cat("_NORM");
  }));

// === Year-to-year change ===
function yearDiff(prevYear, nextYear) {
  var prev = tssNorm.select("TSS_" + prevYear + "_NORM");
  var next = tssNorm.select("TSS_" + nextYear + "_NORM");
  return next.subtract(prev).rename("TSS_DIFF_" + prevYear + "_" + nextYear);
}

var diffs = [];
for (var i = 0; i < years.length - 1; i++) {
  diffs.push(yearDiff(years[i], years[i + 1]));
}
var diffStack = ee.ImageCollection.fromImages(diffs).toBands();

// === Classify turbidity change ===
var highThresh = 0.3;
var modThresh = 0.15;
var latestDiff = yearDiff(years[years.length - 2], years[years.length - 1]);
var highIncrease = latestDiff.gt(highThresh);
var modIncrease = latestDiff.gt(modThresh).and(latestDiff.lte(highThresh));

var changeClass = ee.Image(0)
  .where(modIncrease, 1)
  .where(highIncrease, 2)
  .rename("CHANGE_CLASS");

// === Persistent increase across multiple years ===
var persistent = diffStack.gt(modThresh).reduce(ee.Reducer.sum()).gte(3);

// === Map layers ===
Map.addLayer(annualCol.first().select("TSS_.*"), {min: 0, max: 300}, "TSS proxy (first year)");
Map.addLayer(latestDiff, {min: -0.3, max: 0.3, palette: ["blue", "white", "red"]}, "Latest TSS change");
Map.addLayer(changeClass, {min: 0, max: 2, palette: ["000000", "FFA500", "FF0000"]}, "Change class");
Map.addLayer(persistent.selfMask(), {palette: ["FF00FF"]}, "Persistent increases");

// === Optional spatial association (not causal) ===
// Intersect hotspots with bare soil near rivers as a spatial hint only.
var bareSoil = ee.Image("ESA/WorldCover/v200/2021").eq(60);
var riverMask = ee.Image("JRC/GSW1_4/Occurrence").gte(70);
var nearRiverBare = bareSoil.and(riverMask.focal_max(500));
var assocMask = persistent.and(nearRiverBare);
Map.addLayer(assocMask.selfMask(), {palette: ["00FFFF"]}, "Persistent + near river bare soil");

// === Exports ===
Export.image.toDrive({
  image: tssNorm,
  description: "turbidity_annual_normalized",
  folder: "GEE_exports",
  fileNamePrefix: "turbidity_annual_normalized",
  region: aoi,
  scale: 30,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: diffStack,
  description: "turbidity_yearly_diffs",
  folder: "GEE_exports",
  fileNamePrefix: "turbidity_yearly_diffs",
  region: aoi,
  scale: 30,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: changeClass,
  description: "turbidity_change_class",
  folder: "GEE_exports",
  fileNamePrefix: "turbidity_change_class",
  region: aoi,
  scale: 30,
  maxPixels: 1e13
});

var hotspotVectors = persistent.selfMask().reduceToVectors({
  geometry: aoi,
  scale: 30,
  geometryType: "polygon",
  maxPixels: 1e13
});

Export.table.toDrive({
  collection: hotspotVectors,
  description: "turbidity_hotspots_vectors",
  folder: "GEE_exports",
  fileNamePrefix: "turbidity_hotspots_vectors"
});
