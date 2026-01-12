// Setup
var region = ee.FeatureCollection("FAO/GAUL/2015/level1")
  .filter(ee.Filter.eq("ADM1_NAME", "Cabanas"))
  .filter(ee.Filter.eq("ADM0_NAME", "El Salvador"));

// Landsat 2 imagery
var landsat2 = ee.ImageCollection("LANDSAT/LM02/C02/T1")
  .filterDate("1975-01-01", "1977-12-31")
  .filterBounds(region)
  .median();

// NDVI mask
var ndvi = landsat2.normalizedDifference(["B4", "B2"]).rename("NDVI");
var ndviMask = ndvi.gt(0.3);

Map.centerObject(region, 9);
Map.addLayer(ndvi.updateMask(ndviMask), {min: 0, max: 1}, "NDVI mask");
