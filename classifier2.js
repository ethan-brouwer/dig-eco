// ASM candidate mapping (no ground truth). Results are likelihood-only, not confirmed mining.

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
var slopeMaxDeg = 20;
var riverBufferM = 1000;
var asmThreshold = 0.7;
var useSentinel1 = true;

// === Masks and helpers ===
function maskS2Aggressive(img) {
  var scl = img.select("SCL");
  var qa = img.select("QA60");
  var cloud = qa.bitwiseAnd(1 << 10).neq(0).or(qa.bitwiseAnd(1 << 11).neq(0));
  var shadow = scl.eq(3);
  var cirrus = scl.eq(10);
  var snow = scl.eq(11);
  var cloudMask = cloud.or(shadow).or(cirrus).or(snow);
  var clean = img.updateMask(cloudMask.not());

  // Conservative spectral cleanup (reduces false positives, may miss some ASM)
  var mndwi = clean.normalizedDifference(["B3", "B11"]);
  var bsi = clean.expression(
    "((SWIR + RED) - (NIR + BLUE)) / ((SWIR + RED) + (NIR + BLUE))",
    {SWIR: clean.select("B11"), RED: clean.select("B4"), NIR: clean.select("B8"), BLUE: clean.select("B2")}
  );
  var notTurbidWater = mndwi.lt(0.4);
  var notExtremeBare = bsi.lt(0.6);
  return clean.updateMask(notTurbidWater).updateMask(notExtremeBare);
}

function addS2Indices(img) {
  var ndvi = img.normalizedDifference(["B8", "B4"]).rename("NDVI");
  var ndwi = img.normalizedDifference(["B3", "B8"]).rename("NDWI");
  var mndwi = img.normalizedDifference(["B3", "B11"]).rename("MNDWI");
  var bsi = img.expression(
    "((SWIR + RED) - (NIR + BLUE)) / ((SWIR + RED) + (NIR + BLUE))",
    {SWIR: img.select("B11"), RED: img.select("B4"), NIR: img.select("B8"), BLUE: img.select("B2")}
  ).rename("BSI");
  var ndbi = img.normalizedDifference(["B11", "B8"]).rename("NDBI");
  return img.addBands([ndvi, ndwi, mndwi, bsi, ndbi]);
}

function seasonalComposite(year) {
  var start = ee.Date(year + dryStart);
  var end = ee.Date((year + 1) + dryEnd);
  var s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
    .filterBounds(aoi)
    .filterDate(start, end)
    .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", cloudMax))
    .map(maskS2Aggressive)
    .map(addS2Indices)
    .median()
    .clip(aoi);
  return s2;
}

function addS1Bands(img) {
  var vv = img.select("VV");
  var vh = img.select("VH");
  var ratio = vv.divide(vh).rename("VVVH");
  return img.addBands(ratio);
}

function s1Composite(start, end) {
  var s1 = ee.ImageCollection("COPERNICUS/S1_GRD")
    .filterBounds(aoi)
    .filterDate(start, end)
    .filter(ee.Filter.eq("instrumentMode", "IW"))
    .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VV"))
    .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VH"))
    .select(["VV", "VH"])
    .median()
    .clip(aoi);
  return addS1Bands(s1);
}

// === Spatial constraints ===
var slope = ee.Terrain.slope(ee.Image("USGS/SRTMGL1_003"));
var slopeMask = slope.lt(slopeMaxDeg);

var rivers = ee.FeatureCollection("WWF/HydroSHEDS/v1/HydroRIVERS")
  .filterBounds(aoi);
var riverMask = rivers
  .map(function(f) { return f.buffer(riverBufferM); })
  .union(1)
  .geometry()
  .bounds()
  .intersection(aoi, 1)
  .buffer(0);

var landcover = ee.Image("ESA/WorldCover/v200/2021");
var urbanMask = landcover.eq(50).not(); // exclude built-up

// === Build composites ===
var composites = years.map(function(y) { return seasonalComposite(y); });
var s2Stack = ee.ImageCollection.fromImages(composites).toBands();

var ndviLatest = seasonalComposite(years[years.length - 1]).select("NDVI");
var ndviEarly = seasonalComposite(years[0]).select("NDVI");
var ndviChange = ndviLatest.subtract(ndviEarly).rename("NDVI_CHANGE");

var base = s2Stack.addBands(ndviChange);
if (useSentinel1) {
  var s1 = s1Composite(
    ee.Date(years[years.length - 1] + dryStart),
    ee.Date((years[years.length - 1] + 1) + dryEnd)
  );
  base = base.addBands(s1);
}

// === Constraint mask ===
var analysisMask = slopeMask
  .updateMask(urbanMask)
  .updateMask(ee.Image().paint(riverMask, 1));

var features = base.updateMask(analysisMask);

// === Weak / pseudo labels ===
// Replace these with digitized mine polygons/points for positives.
var asmPositives = ee.FeatureCollection([]);

// Negatives from stable vegetation/agriculture
var vegetation = landcover.eq(10).or(landcover.eq(20));
var stableVeg = vegetation.updateMask(analysisMask).selfMask();
var negPoints = stableVeg.sample({
  region: aoi,
  scale: 30,
  numPixels: 1000,
  seed: 42,
  geometries: true
}).map(function(f) { return f.set("class", 0); });

var posPoints = asmPositives.map(function(f) { return f.set("class", 1); });
var trainingPoints = posPoints.merge(negPoints);

// === Train classifier (probability output) ===
var training = features.sampleRegions({
  collection: trainingPoints,
  properties: ["class"],
  scale: 10,
  tileScale: 4
});

var classifier = ee.Classifier.smileRandomForest({numberOfTrees: 300})
  .setOutputMode("PROBABILITY")
  .train({
    features: training,
    classProperty: "class",
    inputProperties: features.bandNames()
  });

// === Apply model ===
var asmProb = features.classify(classifier).rename("ASM_PROB");
var asmMask = asmProb.gt(asmThreshold).selfMask();

// === Post-process ===
var minPatchPx = 8;
var connected = asmMask.connectedPixelCount(100, true);
var asmFiltered = asmMask.updateMask(connected.gte(minPatchPx)).rename("ASM_MASK");

// === Map layers ===
Map.addLayer(ndviLatest, {min: 0, max: 0.8, palette: ["brown", "yellow", "green"]}, "NDVI (latest dry season)");
Map.addLayer(asmProb, {min: 0, max: 1, palette: ["0020FF", "FFFFFF", "FF0000"]}, "ASM probability");
Map.addLayer(asmFiltered, {palette: ["FF5500"]}, "ASM candidates (filtered)");

// === Exports ===
Export.image.toDrive({
  image: asmProb,
  description: "asm_probability",
  folder: "GEE_exports",
  fileNamePrefix: "asm_probability",
  region: aoi,
  scale: 30,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: asmFiltered,
  description: "asm_candidates_mask",
  folder: "GEE_exports",
  fileNamePrefix: "asm_candidates_mask",
  region: aoi,
  scale: 30,
  maxPixels: 1e13
});

var asmVectors = asmFiltered.reduceToVectors({
  geometry: aoi,
  scale: 30,
  geometryType: "polygon",
  maxPixels: 1e13
});

Export.table.toDrive({
  collection: asmVectors,
  description: "asm_candidates_vectors",
  folder: "GEE_exports",
  fileNamePrefix: "asm_candidates_vectors"
});
