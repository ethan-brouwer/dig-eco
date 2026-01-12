var geometry = /* color: #d63000 */ee.Geometry.MultiPolygon(
        [[[[-88.7931508663418, 13.916793390472376],
           [-88.79343518049738, 13.916665820136247],
           [-88.79356929094813, 13.916462748843754],
           [-88.79346736700556, 13.91637162705198],
           [-88.79321792156718, 13.916301333073784],
           [-88.79303016693613, 13.916470559281379]]],
         [[[-88.82113470118902, 13.952888561922594],
           [-88.82153971475027, 13.95264387334345],
           [-88.82141365092657, 13.952388772207915],
           [-88.82100595515631, 13.952229984623925],
           [-88.82079137843512, 13.952276839987954],
           [-88.82080478948019, 13.952680316339785]]]]),
    imageVisParam = {"opacity":1,"bands":["Mine_65_filtered"],"palette":["ff0000"]};
    var gaul2 = ee.FeatureCollection('FAO/GAUL/2015/level2');
var Cabanas = gaul2.filter(ee.Filter.eq('ADM1_NAME', 'Cabanas'));
var ROI = Cabanas.geometry();
Map.centerObject(ROI, 10);
Map.addLayer(ROI, {});

// Dynamic World masks (shorter period, coarser)
var dw = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')
  .filterBounds(ROI)
  .filterDate('2018-01-01','2024-12-31');

var dwMean = dw.select(['water','built','trees'])
  .mean()
  .reproject({crs: 'EPSG:4326', scale: 30})
  .clip(ROI);

var waterProb = dwMean.select('water');
var builtProb = dwMean.select('built');
var treesProb = dwMean.select('trees');

// CRITICAL FIX #1: Separate masks for training vs. final output
// Training mask: exclude minor urban and water
var nonUrbanNonWaterMask = waterProb.lt(0.25).and(builtProb.lt(0.25));
var forestMask = treesProb.gt(0.8);

// OUTPUT mask: exclude MAJOR urban areas and water bodies
var majorUrbanMask = builtProb.gt(0.6);  // Areas with >50% built probability
var waterMask = waterProb.gt(0.3);       // Areas with >30% water probability
var outputMask = majorUrbanMask.not().and(waterMask.not());

// Visualize masks to verify
Map.addLayer(majorUrbanMask.selfMask(), {palette:['purple']}, 'Major urban areas (masked)', false);
Map.addLayer(waterMask.selfMask(), {palette:['blue']}, 'Water bodies (masked)', false);
Map.addLayer(outputMask.selfMask(), {palette:['gray'], opacity: 0.3}, 'Valid prediction area', false);

function maskS2(img) {
  var scl = img.select('SCL');
  var mask = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
  return img.updateMask(mask);
}

function addIndices(img) {
  var ndvi = img.normalizedDifference(['B8','B4']).rename('NDVI');
  var bsi  = img.expression(
    '((SWIR + RED) - (NIR + BLUE)) / ((SWIR + RED) + (NIR + BLUE))',
    {
      'SWIR': img.select('B11'),
      'RED' : img.select('B4'),
      'NIR' : img.select('B8'),
      'BLUE': img.select('B2')
    }
  ).rename('BSI').clamp(-2, 2);
  var nbi = img.expression(
    '(RED * SWIR) / NIR',
    {'RED': img.select('B4'), 'SWIR': img.select('B11'), 'NIR': img.select('B8')}
  ).rename('NBI');
  return img.addBands([ndvi, bsi, nbi]);
}

function makeDrySeasonComposite(startDate, endDate, suffix) {
  var col = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(ROI)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 80))
    .map(maskS2)
    .map(addIndices);
  
  var count = col.size();
  print('Images in ' + suffix + ':', count);
  
  // ROBUST FIX: Use median and unmask to fill gaps with neutral values
  var comp = col.median()
    .select(['NDVI','BSI','NBI'])
    .unmask(0)  // Fill missing pixels with 0 (neutral value)
    .rename(['NDVI'+suffix, 'BSI'+suffix, 'NBI'+suffix])
    .clip(ROI);
  
  return comp;
}

function harmonizeLandsat(img) {
  var slopes = ee.Image.constant([0.9785, 0.9542, 0.9825, 1.0073, 1.0171, 0.9949]);
  var intercepts = ee.Image.constant([-0.0095, -0.0016, -0.0022, -0.0021, -0.0030, 0.0029]);
  var harmonized = img.select(['B2','B3','B4','B5','B6','B7'])
    .multiply(slopes).add(intercepts)
    .multiply(10000).int16();
  return img.addBands(harmonized, null, true);
}

function addIndicesLandsat(img) {
  var ndvi = img.normalizedDifference(['B5','B4']).rename('NDVI');
  var bsi  = img.expression(
    '((SWIR + RED) - (NIR + BLUE)) / ((SWIR + RED) + (NIR + BLUE))',
    {
      'SWIR': img.select('B6'),
      'RED' : img.select('B4'),
      'NIR' : img.select('B5'),
      'BLUE': img.select('B2')
    }
  ).rename('BSI').clamp(-2, 2);
  var nbi = img.expression(
    '(RED * SWIR) / NIR',
    {'RED': img.select('B4'), 'SWIR': img.select('B6'), 'NIR': img.select('B5')}
  ).rename('NBI');
  return img.addBands([ndvi, bsi, nbi]);
}

function maskCloudsLandsat(img) {
  var qa = img.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 4).eq(0));
  return img.updateMask(mask);
}

function makeLandsatComposite(collection, startDate, endDate, suffix, bandMap) {
  var col = ee.ImageCollection(collection)
    .filterBounds(ROI)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.lt('CLOUD_COVER', 80))
    .map(maskCloudsLandsat)
    .map(function(img) {
      return img.select(bandMap.from, bandMap.to);
    })
    .map(harmonizeLandsat)
    .map(addIndicesLandsat);
  
  var count = col.size();
  print('Landsat images in ' + suffix + ':', count);
  
  // ROBUST FIX: Use median and unmask
  var comp = col.median()
    .select(['NDVI','BSI','NBI'])
    .unmask(0)  // Fill missing pixels with 0
    .rename(['NDVI'+suffix, 'BSI'+suffix, 'NBI'+suffix])
    .clip(ROI);
  return comp;
}

// ALTERNATIVE APPROACH: Create data quality mask to identify valid pixels
var l8map = {from: ['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7'], to: ['B2','B3','B4','B5','B6','B7']};

var dry1718_quality = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(ROI)
  .filterDate('2017-11-01','2018-04-30')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 80))
  .map(maskS2)
  .select('B8')
  .count()
  .gte(1);  // Pixels with at least 1 valid observation

print('Data quality check - pixels with valid 1718 data:', 
      dry1718_quality.reduceRegion({
        reducer: ee.Reducer.sum(),
        geometry: ROI,
        scale: 100,
        maxPixels: 1e9
      }));

// Create overall data quality mask
var dataQualityMask = dry1718_quality;

// Visualize data gaps
Map.addLayer(dataQualityMask.not().selfMask(), 
  {palette:['red'], opacity: 0.5}, 
  'Data gaps (1718 missing)', false);

var dry1314 = makeLandsatComposite('LANDSAT/LC08/C02/T1_L2', '2013-11-01','2014-04-30','_1314', l8map);
var dry1415 = makeLandsatComposite('LANDSAT/LC08/C02/T1_L2', '2014-11-01','2015-04-30','_1415', l8map);
var dry1516 = makeLandsatComposite('LANDSAT/LC08/C02/T1_L2', '2015-11-01','2016-04-30','_1516', l8map);
var dry1617 = makeLandsatComposite('LANDSAT/LC08/C02/T1_L2', '2016-11-01','2017-04-30','_1617', l8map);
var dry1718 = makeDrySeasonComposite('2017-11-01','2018-04-30','_1718');
var dry1819 = makeDrySeasonComposite('2018-11-01','2019-04-30','_1819');
var dry1920 = makeDrySeasonComposite('2019-11-01','2020-04-30','_1920');
var dry2021 = makeDrySeasonComposite('2020-11-01','2021-04-30','_2021');
var dry2122 = makeDrySeasonComposite('2021-11-01','2022-04-30','_2122');
var dry2223 = makeDrySeasonComposite('2022-11-01','2023-04-30','_2223');
var dry2324 = makeDrySeasonComposite('2023-11-01','2024-04-30','_2324');

// OPTION 1: Stack with all time periods (fills gaps with 0)
var stackedImageFull = dry1314
  .addBands([dry1415, dry1516, dry1617, dry1718, dry1819, dry1920, 
             dry2021, dry2122, dry2223, dry2324])
  .clip(ROI);

// OPTION 2: Stack WITHOUT the problematic 1718 period
var stackedImageRobust = dry1314
  .addBands([dry1415, dry1516, dry1617, dry1819, dry1920,  // Skip dry1718
             dry2021, dry2122, dry2223, dry2324])
  .clip(ROI);

print('Full stack bands (33):', stackedImageFull.bandNames().length());
print('Robust stack bands (30 - no 1718):', stackedImageRobust.bandNames().length());

// USE THE ROBUST STACK for classification
var stackedImage = stackedImageRobust;

// Training stack WITH mask
var stackedTrain = stackedImage.updateMask(nonUrbanNonWaterMask);

Map.addLayer(
  stackedImage.select(['NDVI_2223','BSI_2223','NBI_2223']),
  {min:[0,0,0], max:[0.8,0.5,2]},
  'Recent indices',
  false
);

// Create a coverage mask showing valid data areas
var coverageMask = stackedImage.select(0).mask();
Map.addLayer(coverageMask.not().selfMask(), 
  {palette:['red'], opacity: 0.4}, 
  'No data areas (masked out)', false);

// Visualize recent NDVI to check coverage
Map.addLayer(dry2223.select('NDVI_2223'),
  {min:0, max:0.8, palette:['brown','yellow','green']},
  'NDVI 2022-23 (check coverage)', false);

var mineSites = ee.FeatureCollection('projects/metalminingpersonalcopy/assets/MStraining')
  .filterBounds(ROI)
  .map(function(f){ return f.buffer(30).set('class', 1); });

var minePoints = mineSites;

print('Number of mine training points:', minePoints.size());

var forestSample = forestMask.selfMask().clip(ROI);
var forestPoints = forestSample.sample({
  region: ROI,
  scale: 30,
  numPixels: 500,
  seed: 44,
  geometries: true
}).map(function(f) {
  return ee.Feature(f.geometry(), {'class': 0});
});

var nonMineSample = nonUrbanNonWaterMask.selfMask()
  .updateMask(forestMask.not())
  .clip(ROI);
  
var nonMinePoints = nonMineSample.sample({
  region: ROI,
  scale: 30,
  numPixels: 1000,
  seed: 42,
  geometries: true
}).map(function(f) {
  return ee.Feature(f.geometry(), {'class': 0});
});

// STRATIFIED TRAIN / VALIDATION SPLIT (replace your old block with this)

// Merge all non-mine points into one collection
var nonMineAll = forestPoints.merge(nonMinePoints);

// Split mine points
var mineSplit = minePoints.randomColumn('rand_mine', 43);
var mineTrain = mineSplit.filter('rand_mine <= 0.7');
var mineVal   = mineSplit.filter('rand_mine > 0.7');

// Split non-mine points
var nonMineSplit = nonMineAll.randomColumn('rand_nonmine', 44);
var nonMineTrain = nonMineSplit.filter('rand_nonmine <= 0.7');
var nonMineVal   = nonMineSplit.filter('rand_nonmine > 0.7');

// Merge into final train / validation sets
var trainPoints = mineTrain.merge(nonMineTrain);
var valPoints   = mineVal.merge(nonMineVal);

print('Training set - mines:', mineTrain.size());
print('Training set - non-mines:', nonMineTrain.size());
print('Validation set - mines:', mineVal.size());
print('Validation set - non-mines:', nonMineVal.size());

print('Training class distribution:', trainPoints.aggregate_histogram('class'));
print('Validation class distribution:', valPoints.aggregate_histogram('class'));

// Sample training data from stackedTrain
var training = stackedTrain.sampleRegions({
  collection: trainPoints,
  properties: ['class'],
  scale: 10,
  tileScale: 4
});


print('Training samples extracted:', training.size());

var posTrain = training.filter('class == 1');
var negTrain = training.filter('class == 0');
var posSize = posTrain.size();

print('Positive training samples:', posSize);
print('Negative training samples (before balancing):', negTrain.size());

var negLimit = posSize.multiply(2.5).toInt();
negTrain = negTrain.limit(negLimit);

var trainBalanced = posTrain.merge(negTrain);
print('Balanced training size:', trainBalanced.size());

var featureBands = stackedTrain.bandNames();
print('Total features:', featureBands.length());

var classifier = ee.Classifier.smileRandomForest({numberOfTrees: 500})
  .train({
    features: trainBalanced,
    classProperty: 'class',
    inputProperties: featureBands
  });

// Classify the UNMASKED stackedImage
var classified = stackedImage.classify(classifier);

// CRITICAL FIX #5: Apply output mask to remove urban/water false positives
var classifiedMasked = classified.updateMask(outputMask);

var trainConf = classifier.confusionMatrix();
print('Training accuracy:', trainConf.accuracy());
print('Training confusion matrix:', trainConf);

var validation = stackedTrain.sampleRegions({
  collection: valPoints,
  properties: ['class'],
  scale: 10,
  tileScale: 4
});

var validated = validation.classify(classifier);
var valConf = validated.errorMatrix('class','classification');
print('Validation accuracy:', valConf.accuracy());
print('Validation confusion matrix:', valConf);
print('Validation producers accuracy:', valConf.producersAccuracy());
print('Validation consumers accuracy:', valConf.consumersAccuracy());

var importance = classifier.explain();
var importanceDict = ee.Dictionary(importance.get('importance'));
var importanceList = importanceDict.keys().zip(importanceDict.values());
var sortedImportance = ee.List(importanceList)
  .sort(ee.List(importanceList).map(function(x){ return ee.List(x).get(1); }))
  .reverse();
print('Top 15 features:', sortedImportance.slice(0, 15));

// Show both raw and masked predictions
var minesOnly = classified.eq(1).selfMask();
var minesOnlyMasked = classifiedMasked.eq(1).selfMask();

Map.addLayer(minesOnly, {palette:['red'], opacity: 0.5}, 'Raw predictions (with false positives)', false);
Map.addLayer(minesOnlyMasked, {palette:['FF0000']}, 'Predicted mines (urban/water masked)', true);

// Probability classifier
var probClassifier = ee.Classifier.smileRandomForest({numberOfTrees: 500})
  .setOutputMode('PROBABILITY')
  .train({
    features: trainBalanced,
    classProperty: 'class',
    inputProperties: featureBands
  });

var mineProb = stackedImage.classify(probClassifier).rename('mineProb');
var mineProbMasked = mineProb.updateMask(outputMask);

Map.addLayer(mineProbMasked, {min:0, max:1, palette:['0000FF','FFFFFF','FF0000']}, 'Mine probability (masked)', true);

// Check probability stats
var probStats = mineProb.reduceRegion({
  reducer: ee.Reducer.minMax(),
  geometry: ROI,
  scale: 30,
  maxPixels: 1e9,
  tileScale: 4
});
print('Probability range:', probStats);

// Area calculations using MASKED predictions
var mineAreaImage = classifiedMasked.eq(1).selfMask().multiply(ee.Image.pixelArea());
var mineAreaStats = mineAreaImage.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: ROI,
  scale: 30,
  maxPixels: 1e10,
  bestEffort: true,
  tileScale: 4
});

print('Estimated mine area - MASKED (hectares):',
      ee.Number(mineAreaStats.get('classification')).divide(10000));

var mine50 = mineProbMasked.gt(0.5).selfMask();
var mine70 = mineProbMasked.gt(0.7).selfMask();
var mine80 = mineProbMasked.gt(0.8).selfMask();


Map.addLayer(mine50, {palette:['FFA500']}, 'Mine prob > 0.5 (masked)', false);
Map.addLayer(mine70, {palette:['FFA500']}, 'Mine prob > 0.7 (masked)', false);
Map.addLayer(mine80, {palette:['00FF00']}, 'Mine prob > 0.8 (masked)', false);

var mine50Area = mine50.multiply(ee.Image.pixelArea()).reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: ROI,
  scale: 30,
  maxPixels: 1e10,
  bestEffort: true,
  tileScale: 4
});

var mine80Area = mine80.multiply(ee.Image.pixelArea()).reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: ROI,
  scale: 30,
  maxPixels: 1e10,
  bestEffort: true,
  tileScale: 4
});

var mine70Area = mine70.multiply(ee.Image.pixelArea()).reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: ROI,
  scale: 30,
  maxPixels: 1e10,
  bestEffort: true,
  tileScale: 4
});
print('Mine area prob>0.5 - MASKED (ha):',
      ee.Number(mine80Area.get('mineProb')).divide(10000));
      
      
      
// === APPLY NON-URBAN / NON-WATER MASK ===
var mineProbMasked = mineProb.updateMask(nonUrbanNonWaterMask);

// === 65% threshold ===
var mineThresh65 = mineProbMasked.gt(0.65).selfMask().rename('Mine_65');

// === Remove tiny isolated patches ===
var minSize = 5;  
var connected = mineThresh65.connectedPixelCount(100, true);
var mineFiltered = mineThresh65.updateMask(connected.gte(minSize)).rename('Mine_65_filtered');

// === Map preview ===
Map.addLayer(mineFiltered, {palette:['ff0000']}, 'Mines > 65% masked', true);

// Convert mine training points to small polygons (control size here)
var minePoly = minePoints.map(function(f) {
  return ee.Feature(f.buffer(15)).set('class', 1);   // 15 m buffer → small polygon
});

// Convert non-mine points to polygons (negative samples)
var nonMinePoly = forestPoints.merge(nonMinePoints).map(function(f) {
  return ee.Feature(f.buffer(15)).set('class', 0);
});

// Add to map
Map.addLayer(minePoly, {color: 'red'}, 'Mine Polygons (class=1)');
Map.addLayer(nonMinePoly, {color: 'blue'}, 'Non-Mine Polygons (class=0)');


Export.image.toDrive({
  image: mine80,
  description: 'mine_prob_gt_08_masked',
  folder: 'GEE_exports',
  fileNamePrefix: 'mine_prob_gt_08_masked',
  region: ROI,
  scale: 30,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});


Export.image.toDrive({
  image: mineFiltered,
  description: 'mine_prob_gt_065_masked',
  folder: 'GEE_exports',
  fileNamePrefix: 'mine_prob_gt_065_masked',
  region: ROI,
  scale: 30,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

