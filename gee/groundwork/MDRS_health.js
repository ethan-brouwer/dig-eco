/*
  FILE: gee/groundwork/MDRS_health.js
  PURPOSE: Monitor dry-season index behavior within 1 km of MRDS sites.
  INPUTS: Landsat 8/9 L2, MRDS CSV asset, GAUL boundaries, WorldCover, JRC water.
  OUTPUTS: Per-site normalized index charts and optional clipped map layers.
  AUTHOR: Ethan Brouwer
  LAST MODIFIED: 2026-01-16
*/

var L8_ID = "LANDSAT/LC08/C02/T1_L2";
var L9_ID = "LANDSAT/LC09/C02/T1_L2";
var GAUL0_ID = "FAO/GAUL/2015/level0";
var GAUL1_ID = "FAO/GAUL/2015/level1";
var WC_ID = "ESA/WorldCover/v200/2021";
var JRC_ID = "JRC/GSW1_4/Occurrence";

var mrdsTable = ee.FeatureCollection("projects/metalminingpersonalcopy/assets/USGS_MRDS");
var nameField = "Name";
var descField = "description";
var xField = "X";
var yField = "Y";

var startYear = 2015;
var endYear = ee.Date(Date.now()).get("year");
var dryMonths = [11, 12, 1, 2, 3, 4];
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

function scaleL2(img) {
  var scaled = img.select([
    bands.blue, bands.green, bands.red, bands.nir, bands.swir1, bands.swir2
  ]).multiply(reflMult).add(reflAdd);
  return img.addBands(scaled, null, true);
}

function maskLandsatL2(img) {
  var qa = img.select("QA_PIXEL");
  var cloud = qa.bitwiseAnd(1 << 3).neq(0);
  var shadow = qa.bitwiseAnd(1 << 4).neq(0);
  return img.updateMask(cloud.or(shadow).not());
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

function extractCommodity(desc) {
  var s = ee.String(desc);
  var m = s.match("commod1</th><td>([^<]*)");
  return ee.Algorithms.If(m.size().gt(1), ee.String(m.get(1)), "unknown");
}

var mrds = mrdsTable.map(function (f) {
  var x = ee.Number(f.get(xField));
  var y = ee.Number(f.get(yField));
  var geom = ee.Geometry.Point([x, y]);
  var commod = extractCommodity(f.get(descField));
  return ee.Feature(geom, f.toDictionary()).set("commod1", commod);
});

var bufferMeters = 1000;
var buffers = mrds.map(function (f) {
  return f.buffer(bufferMeters).copyProperties(f).set("buffer_m", bufferMeters);
});

Map.addLayer(buffers.style({color: "FFFFFF", width: 2, fillColor: "00000000"}), {}, "MRDS buffers (1 km)", false);
Map.addLayer(mrds.style({color: "FFAA00", pointSize: 6}), {}, "MRDS sites", false);

var worldcover = ee.Image(WC_ID);
var urbanMask = worldcover.neq(50);
var waterOcc = ee.Image(JRC_ID).gte(70);
var waterMask = waterOcc.not();
var analysisMask = urbanMask.and(waterMask);

var l8 = ee.ImageCollection(L8_ID);
var l9 = ee.ImageCollection(L9_ID);

var baseCollection = l8.merge(l9)
  .filterBounds(elsal)
  .filter(ee.Filter.lt("CLOUD_COVER", cloudMax))
  .map(scaleL2)
  .map(maskLandsatL2)
  .map(addIndices)
  .map(function (img) {
    return img.updateMask(analysisMask);
  });

var indexBands = [
  "NDVI", "GNDVI", "NDWI", "MNDWI", "NDMI", "NDBI", "NDTI",
  "BSI", "SAVI", "IOI", "KAOLINITE", "CLAY", "FERROUS"
];

function monthComposite(year, month) {
  var start = ee.Date.fromYMD(year, month, 1);
  var end = start.advance(1, "month");
  return baseCollection
    .filterDate(start, end)
    .median()
    .set({
      year: year,
      month: month,
      date: start.format("YYYY-MM")
    });
}

function buildMonthlyCollection() {
  var years = ee.List.sequence(startYear, endYear);
  var months = ee.List(dryMonths);
  var images = years.map(function (y) {
    y = ee.Number(y);
    return months.map(function (m) {
      return monthComposite(y, m);
    });
  }).flatten();
  return ee.ImageCollection.fromImages(images);
}

var monthly = buildMonthlyCollection();

function buildSeries(geom) {
  return ee.FeatureCollection(monthly.map(function (img) {
    var stats = img.select(indexBands).reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: geom,
      scale: 30,
      bestEffort: true,
      maxPixels: 1e9
    });
    return ee.Feature(null, stats).set({
      date: img.get("date")
    });
  }));
}

function normalizeSeries(fc) {
  var mins = ee.Dictionary(indexBands.map(function (b) {
    b = ee.String(b);
    return ee.List([b, fc.aggregate_min(b)]);
  }).flatten());

  var maxs = ee.Dictionary(indexBands.map(function (b) {
    b = ee.String(b);
    return ee.List([b, fc.aggregate_max(b)]);
  }).flatten());

  return fc.map(function (f) {
    var props = ee.Dictionary(indexBands.map(function (b) {
      b = ee.String(b);
      var val = ee.Number(f.get(b));
      var min = ee.Number(mins.get(b));
      var max = ee.Number(maxs.get(b));
      var denom = max.subtract(min);
      var norm = ee.Algorithms.If(denom.neq(0), val.subtract(min).divide(denom), 0);
      return ee.List([b, norm]);
    }).flatten());
    return f.set(props);
  });
}

var chartPanel = ui.Panel({style: {position: "bottom-right", width: "420px"}});
var selector = ui.Select({
  items: mrds.aggregate_array(nameField).distinct(),
  placeholder: "Select site",
  onChange: function (name) {
    var site = mrds.filter(ee.Filter.eq(nameField, name)).first();
    var geom = ee.Feature(site).geometry().buffer(bufferMeters);

    var series = normalizeSeries(buildSeries(geom));

    var chart = ui.Chart.feature.byFeature({
      features: series,
      xProperty: "date",
      yProperties: indexBands
    }).setOptions({
      title: name + " (dry season monthly medians, normalized)",
      lineWidth: 1,
      pointSize: 2,
      vAxis: {title: "Normalized value"},
      hAxis: {title: "Year-Month"}
    });

    chartPanel.clear();
    chartPanel.add(ui.Label("Site: " + name));
    chartPanel.add(ui.Label("Dry season: Nov–Apr (Feb often most stable)"));
    chartPanel.add(chart);
  }
});

chartPanel.add(ui.Label("MRDS site selector"));
chartPanel.add(selector);
Map.add(chartPanel);

// Optional clipped layers for a selected buffer
// var exampleSite = ee.Feature(mrds.first()).geometry().buffer(bufferMeters);
// Map.addLayer(monthly.first().select("NDVI").clip(exampleSite), {min: 0, max: 0.8, palette: ["brown", "yellow", "green"]}, "NDVI (buffer)");
