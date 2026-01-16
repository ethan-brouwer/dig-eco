/*
  FILE: gee/groundwork/GEEstarter.js
  PURPOSE: Baseline visualization for El Salvador with Landsat indices.
  INPUTS: Landsat 8/9 L2, GAUL boundaries, user-defined date/month filters.
  OUTPUTS: Map-ready composites and index layers (commented out by default).
  AUTHOR: Ethan Brouwer
  LAST MODIFIED: 2026-01-16
*/

var L8_ID = "LANDSAT/LC08/C02/T1_L2";
var L9_ID = "LANDSAT/LC09/C02/T1_L2";
var GAUL0_ID = "FAO/GAUL/2015/level0";
var GAUL1_ID = "FAO/GAUL/2015/level1";

var startDate = "2019-01-01";
var endDate = "2024-12-31";
var months = [1, 2, 3, 4];
var cloudMax = 60;

var reflMult = 0.0000275;
var reflAdd = -0.2;
var saviL = 0.5;

var bands = {
  blue: "SR_B2",
  green: "SR_B3",
  red: "SR_B4",
  nir: "SR_B5",
  swir1: "SR_B6",
  swir2: "SR_B7"
};

var elsal = ee.FeatureCollection(GAUL0_ID)
  .filter(ee.Filter.eq("ADM0_NAME", "El Salvador"));

var departments = ee.FeatureCollection(GAUL1_ID)
  .filter(ee.Filter.eq("ADM0_NAME", "El Salvador"));

Map.setOptions("SATELLITE");
Map.centerObject(elsal, 8);
// Map.addLayer(elsal.style({color: "FFFFFF", width: 2, fillColor: "00000000"}), {}, "El Salvador");
// Map.addLayer(departments.style({color: "00FFFF", width: 1, fillColor: "00000000"}), {}, "Departments");

function maskLandsatL2(img) {
  var qa = img.select("QA_PIXEL");
  var cloud = qa.bitwiseAnd(1 << 3).neq(0);
  var shadow = qa.bitwiseAnd(1 << 4).neq(0);
  var mask = cloud.or(shadow).not();
  return img.updateMask(mask);
}

function scaleL2(img) {
  var scaled = img.select([
    bands.blue, bands.green, bands.red, bands.nir, bands.swir1, bands.swir2
  ]).multiply(reflMult).add(reflAdd);
  return img.addBands(scaled, null, true);
}

function addIndices(img) {
  var blue = img.select(bands.blue);
  var green = img.select(bands.green);
  var red = img.select(bands.red);
  var nir = img.select(bands.nir);
  var swir1 = img.select(bands.swir1);
  var swir2 = img.select(bands.swir2);

  var ndvi = img.normalizedDifference([bands.nir, bands.red]).rename("NDVI");
  var gndvi = img.normalizedDifference([bands.nir, bands.green]).rename("GNDVI");
  var ndwi = img.normalizedDifference([bands.green, bands.nir]).rename("NDWI");
  var mndwi = img.normalizedDifference([bands.green, bands.swir1]).rename("MNDWI");
  var ndmi = img.normalizedDifference([bands.nir, bands.swir1]).rename("NDMI");
  var ndbi = img.normalizedDifference([bands.swir1, bands.nir]).rename("NDBI");
  var ndti = img.normalizedDifference([bands.red, bands.green]).rename("NDTI");

  var bsi = img.expression(
    "((SWIR + RED) - (NIR + BLUE)) / ((SWIR + RED) + (NIR + BLUE))",
    {SWIR: swir1, RED: red, NIR: nir, BLUE: blue}
  ).rename("BSI");

  var savi = img.expression(
    "((NIR - RED) / (NIR + RED + L)) * (1 + L)",
    {NIR: nir, RED: red, L: saviL}
  ).rename("SAVI");

  var ioi = red.divide(blue).rename("IOI");
  var kaolinite = swir1.divide(swir2).rename("KAOLINITE");
  var clay = swir1.divide(swir2).rename("CLAY");
  var ferrous = swir1.divide(nir).rename("FERROUS");

  return img.addBands([
    ndvi, gndvi, ndwi, mndwi, ndmi, ndbi, ndti, bsi, savi, ioi, kaolinite, clay, ferrous
  ]);
}

function addMonth(img) {
  var m = ee.Date(img.get("system:time_start")).get("month");
  return img.set("month", m);
}

var l8 = ee.ImageCollection(L8_ID);
var l9 = ee.ImageCollection(L9_ID);

var collection = l8.merge(l9)
  .filterBounds(elsal)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.lt("CLOUD_COVER", cloudMax))
  .map(scaleL2)
  .map(maskLandsatL2)
  .map(addIndices)
  .map(addMonth)
  .filter(ee.Filter.inList("month", months));

var composite = collection.median().clip(elsal);

var ndvi = composite.select("NDVI");
var gndvi = composite.select("GNDVI");
var ndwi = composite.select("NDWI");
var mndwi = composite.select("MNDWI");
var ndmi = composite.select("NDMI");
var ndbi = composite.select("NDBI");
var ndti = composite.select("NDTI");
var bsi = composite.select("BSI");
var savi = composite.select("SAVI");
var ioi = composite.select("IOI");
var kaolinite = composite.select("KAOLINITE");
var clay = composite.select("CLAY");
var ferrous = composite.select("FERROUS");

// Map.addLayer(ndvi, {min: 0, max: 0.8, palette: ["brown", "yellow", "green"]}, "NDVI");
// Map.addLayer(gndvi, {min: 0, max: 0.8, palette: ["brown", "yellow", "green"]}, "GNDVI");
// Map.addLayer(ndwi, {min: -0.5, max: 0.5, palette: ["brown", "white", "blue"]}, "NDWI");
// Map.addLayer(mndwi, {min: -0.5, max: 0.5, palette: ["brown", "white", "blue"]}, "MNDWI");
// Map.addLayer(ndmi, {min: -0.5, max: 0.5, palette: ["brown", "white", "blue"]}, "NDMI");
// Map.addLayer(ndbi, {min: -0.5, max: 0.5, palette: ["white", "gray", "black"]}, "NDBI");
// Map.addLayer(ndti, {min: -0.5, max: 0.5, palette: ["blue", "white", "red"]}, "NDTI");
// Map.addLayer(bsi, {min: -0.5, max: 0.5, palette: ["0000FF", "FFFFFF", "FF0000"]}, "BSI");
// Map.addLayer(savi, {min: 0, max: 0.8, palette: ["brown", "yellow", "green"]}, "SAVI");
// Map.addLayer(ioi, {min: 0.6, max: 2.0, palette: ["blue", "white", "red"]}, "IOI");
// Map.addLayer(kaolinite, {min: 0.6, max: 1.6, palette: ["purple", "white", "orange"]}, "Kaolinite");
// Map.addLayer(clay, {min: 0.6, max: 1.6, palette: ["purple", "white", "orange"]}, "Clay");
// Map.addLayer(ferrous, {min: 0.6, max: 1.6, palette: ["gray", "white", "red"]}, "Ferrous");

// Export.image.toDrive({
//   image: composite.select(["NDVI", "BSI", "NDBI", "IOI", "KAOLINITE"]),
//   description: "elsal_baseline_indices",
//   folder: "GEE_exports",
//   fileNamePrefix: "elsal_baseline_indices",
//   region: elsal.geometry(),
//   scale: 30,
//   crs: "EPSG:4326",
//   maxPixels: 1e13
// });
