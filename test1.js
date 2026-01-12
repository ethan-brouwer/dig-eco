// Setup
var gaul1 = ee.FeatureCollection('FAO/GAUL/2015/level1');
var cOnly = gaul1
  .filter(ee.Filter.eq('ADM0_NAME', 'El Salvador'))
  .filter(ee.Filter.eq('ADM1_NAME', 'Cabanas'));
var roi = cOnly.geometry();
Map.centerObject(roi, 10);

// Latest 4 months of Sentinel-2 SR
var end = ee.Date(Date.now());
var start = end.advance(-4, 'month');

function maskS2(img) {
  var scl = img.select('SCL');
  var mask = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
  return img.updateMask(mask);
}

var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(roi)
  .filterDate(start, end)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
  .map(maskS2)
  .median()
  .clip(roi);

// NDVI layer
var ndvi = s2.normalizedDifference(['B8', 'B4']).rename('NDVI');
Map.addLayer(ndvi, {min: 0, max: 0.8, palette: ['brown', 'yellow', 'green']}, 'NDVI (last 4 months)');
Map.addLayer(roi, {color: 'white'}, 'Cabanas boundary', false);
