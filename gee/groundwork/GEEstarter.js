/*
  FILE: gee/groundwork/GEEstarter.js
  PURPOSE: Baseline visualization for El Salvador with Landsat indices.
  INPUTS: Landsat 8/9 L2, GAUL boundaries, user-defined date/month filters.
  OUTPUTS: Map-ready composites and index layers (commented out by default).
  AUTHOR: Ethan Brouwer
  LAST MODIFIED: 2026-01-16
*/

var startDate = "2019-01-01";
var endDate = "2024-12-31";
var months = [1, 2, 3, 4];
var cloudMax = 60;

var elsal = ee.FeatureCollection("FAO/GAUL/2015/level0")
  .filter(ee.Filter.eq("ADM0_NAME", "El Salvador"));

var departments = ee.FeatureCollection("FAO/GAUL/2015/level1")
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

function addIndices(img) {
  var ndvi = img.normalizedDifference(["SR_B5", "SR_B4"]).rename("NDVI");
  var bsi = img.expression(
    "((SWIR + RED) - (NIR + BLUE)) / ((SWIR + RED) + (NIR + BLUE))",
    {
      SWIR: img.select("SR_B6"),
      RED: img.select("SR_B4"),
      NIR: img.select("SR_B5"),
      BLUE: img.select("SR_B2")
    }
  ).rename("BSI");
  var ndwi = img.normalizedDifference(["SR_B3", "SR_B5"]).rename("NDWI");
  var ndbi = img.normalizedDifference(["SR_B6", "SR_B5"]).rename("NDBI");
  var ioi = img.select("SR_B4").divide(img.select("SR_B3")).rename("IOI");
  var kaolinite = img.select("SR_B6").divide(img.select("SR_B7")).rename("KAOLINITE");
  var gndvi = img.normalizedDifference(["SR_B5", "SR_B3"]).rename("GNDVI");
  return img.addBands([ndvi, bsi, ndwi, ndbi, ioi, kaolinite, gndvi]);
}

function monthFilter(img) {
  var m = ee.Date(img.get("system:time_start")).get("month");
  return ee.List(months).contains(m);
}

var l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2");
var l9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2");

var collection = l8.merge(l9)
  .filterBounds(elsal)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.lt("CLOUD_COVER", cloudMax))
  .filter(monthFilter)
  .map(maskLandsatL2)
  .map(addIndices);

var composite = collection.median().clip(elsal);

var ndvi = composite.select("NDVI");
var bsi = composite.select("BSI");
var ndwi = composite.select("NDWI");
var ndbi = composite.select("NDBI");
var ioi = composite.select("IOI");
var kaolinite = composite.select("KAOLINITE");
var gndvi = composite.select("GNDVI");

// Map.addLayer(ndvi, {min: 0, max: 0.8, palette: ["brown", "yellow", "green"]}, "NDVI");
// Map.addLayer(bsi, {min: -0.5, max: 0.5, palette: ["0000FF", "FFFFFF", "FF0000"]}, "BSI");
// Map.addLayer(ndwi, {min: -0.5, max: 0.5, palette: ["brown", "white", "blue"]}, "NDWI");
// Map.addLayer(ndbi, {min: -0.5, max: 0.5, palette: ["white", "gray", "black"]}, "NDBI");
// Map.addLayer(ioi, {min: 0.6, max: 1.6, palette: ["blue", "white", "red"]}, "IOI");
// Map.addLayer(kaolinite, {min: 0.6, max: 1.4, palette: ["purple", "white", "orange"]}, "Kaolinite");
// Map.addLayer(gndvi, {min: 0, max: 0.8, palette: ["brown", "yellow", "green"]}, "GNDVI");
