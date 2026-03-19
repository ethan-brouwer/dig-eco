/*
  FILE: Phase I/gee_scripts/turbidity/san_sebastian_downstream_polygon_valid_pixels_feb2023.js
  PURPOSE: Quick feasibility check for valid Sentinel-2 pixels inside user-drawn
           downstream polygons for February 2023.

  GEE IMPORTS REQUIRED
  - impact_pool
  - Poly500m
  - Poly1000m
  - Poly1500m
  - Poly2000m
  - Poly2500m
  - Poly3000m
  - Poly3500m
  - Poly4000m
  - Poly4500m
  - Poly5000m
*/

// ================= USER SETTINGS =================
var cloudMax = 60;
var analysisScaleMeters = 10;
var monthStart = ee.Date("2023-02-01");
var monthEnd = monthStart.advance(1, "month");

// ================= HELPERS =================
function toGeometry(obj) {
  var t = ee.String(ee.Algorithms.ObjectType(obj));
  return ee.Geometry(ee.Algorithms.If(
    t.compareTo("FeatureCollection").eq(0),
    ee.FeatureCollection(obj).geometry(),
    ee.Algorithms.If(
      t.compareTo("Feature").eq(0),
      ee.Feature(obj).geometry(),
      ee.Geometry(obj)
    )
  ));
}

function maskS2(img) {
  var scl = img.select("SCL");
  var qa = img.select("QA60");

  var cloud = qa.bitwiseAnd(1 << 10).neq(0)
    .or(qa.bitwiseAnd(1 << 11).neq(0));
  var shadow = scl.eq(3);
  var cirrus = scl.eq(10);
  var snow = scl.eq(11);
  var saturated = scl.eq(1);

  var bad = cloud.or(shadow).or(cirrus).or(snow).or(saturated);
  return img.updateMask(bad.not());
}

function scaleAndSelectS2(img) {
  return img.select(["B2", "B3", "B4", "B8", "QA60", "SCL"], [
    "B2", "B3", "B4", "B8", "QA60", "SCL"
  ]).addBands(
    img.select(["B2", "B3", "B4", "B8"])
      .multiply(0.0001),
    null,
    true
  );
}

function addIndices(img) {
  var ndssi = img.expression(
    "(BLUE - NIR) / (BLUE + NIR)",
    {
      BLUE: img.select("B2"),
      NIR: img.select("B8")
    }
  ).rename("NDSSI");

  var egri = img.expression(
    "GREEN / RED",
    {
      GREEN: img.select("B3"),
      RED: img.select("B4").max(0.0001)
    }
  ).rename("EGRI");

  return img.addBands([ndssi, egri]);
}

function summarizePolygon(feature) {
  feature = ee.Feature(feature);
  var geom = feature.geometry();
  var areaSqm = geom.area(1);

  var validMask = monthImage.select("B4").mask().clip(geom);
  var validPx = ee.Number(ee.Algorithms.If(
    imageCount.gt(0),
    validMask.reduceRegion({
      reducer: ee.Reducer.count(),
      geometry: geom,
      scale: analysisScaleMeters,
      bestEffort: true,
      maxPixels: 1e8
    }).get("B4"),
    0
  ));

  var stats = ee.Dictionary(ee.Algorithms.If(
    validPx.gt(0),
    monthImage.select(["NDSSI", "EGRI"]).updateMask(validMask).reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: geom,
      scale: analysisScaleMeters,
      bestEffort: true,
      maxPixels: 1e8
    }),
    ee.Dictionary({})
  ));

  return feature.set({
    polygon_area_sqm: areaSqm,
    image_count: imageCount,
    valid_px: validPx,
    has_valid_pixels: validPx.gt(0),
    ndssi_mean: stats.get("NDSSI"),
    egri_mean: stats.get("EGRI")
  });
}

// ================= INPUTS =================
var polygons = ee.FeatureCollection([
  ee.Feature(toGeometry(impact_pool), {polygon_id: "impact_pool", distance_m: 0}),
  ee.Feature(toGeometry(Poly500m), {polygon_id: "Poly500m", distance_m: 500}),
  ee.Feature(toGeometry(Poly1000m), {polygon_id: "Poly1000m", distance_m: 1000}),
  ee.Feature(toGeometry(Poly1500m), {polygon_id: "Poly1500m", distance_m: 1500}),
  ee.Feature(toGeometry(Poly2000m), {polygon_id: "Poly2000m", distance_m: 2000}),
  ee.Feature(toGeometry(Poly2500m), {polygon_id: "Poly2500m", distance_m: 2500}),
  ee.Feature(toGeometry(Poly3000m), {polygon_id: "Poly3000m", distance_m: 3000}),
  ee.Feature(toGeometry(Poly3500m), {polygon_id: "Poly3500m", distance_m: 3500}),
  ee.Feature(toGeometry(Poly4000m), {polygon_id: "Poly4000m", distance_m: 4000}),
  ee.Feature(toGeometry(Poly4500m), {polygon_id: "Poly4500m", distance_m: 4500}),
  ee.Feature(toGeometry(Poly5000m), {polygon_id: "Poly5000m", distance_m: 5000})
]);

var aoi = polygons.geometry().bounds().buffer(250);

// ================= COLLECTION =================
var s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  .filterBounds(aoi)
  .filterDate(monthStart, monthEnd)
  .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", cloudMax))
  .map(scaleAndSelectS2)
  .map(maskS2)
  .map(addIndices);

var imageCount = s2.size();
var monthImage = ee.Image(ee.Algorithms.If(
  imageCount.gt(0),
  s2.median().clip(aoi),
  ee.Image.constant([0, 0, 0, 0, 0, 0]).rename(["B2", "B3", "B4", "B8", "NDSSI", "EGRI"]).clip(aoi).selfMask()
));

// ================= OUTPUTS =================
Map.setOptions("SATELLITE");
Map.centerObject(aoi, 15);
Map.addLayer(polygons, {color: "FFB300"}, "Downstream polygons", true);

Map.addLayer(
  monthImage.select(["B4", "B3", "B2"]),
  {min: 0.02, max: 0.3},
  "February 2023 RGB",
  true
);

Map.addLayer(
  monthImage.select("NDSSI"),
  {min: -0.4, max: 0.4, palette: ["#7B3294", "#F7F7F7", "#008837"]},
  "February 2023 NDSSI",
  false
);

var polygonStats = polygons.map(summarizePolygon).sort("distance_m");

print("Month window", monthStart.format("YYYY-MM-dd"), monthEnd.format("YYYY-MM-dd"));
print("Sentinel-2 image count after filtering", imageCount);
print("Valid pixel feasibility by polygon", polygonStats);
print("Polygons with valid pixels", polygonStats.filter(ee.Filter.gt("valid_px", 0)));
print("Polygons with zero valid pixels", polygonStats.filter(ee.Filter.eq("valid_px", 0)));

var validSignalStats = polygonStats.filter(ee.Filter.gt("valid_px", 0));
print("Polygon signal summary for comparison", validSignalStats);

var validPxChart = ui.Chart.feature.byFeature(
  polygonStats,
  "polygon_id",
  ["valid_px"]
).setChartType("ColumnChart").setOptions({
  title: "February 2023 Valid Pixels by Polygon",
  hAxis: {title: "Polygon"},
  vAxis: {title: "Valid pixels at 10 m"},
  colors: ["#FB8C00"]
});
print(validPxChart);

var egriDistanceChart = ui.Chart.feature.byFeature(
  validSignalStats,
  "distance_m",
  ["egri_mean"]
).setChartType("LineChart").setOptions({
  title: "February 2023 EGRI by Distance Downstream",
  hAxis: {title: "Distance from impact pool (m)"},
  vAxis: {title: "Mean EGRI"},
  lineWidth: 2,
  pointSize: 5,
  colors: ["#2E7D32"]
});
print(egriDistanceChart);

var ndssiDistanceChart = ui.Chart.feature.byFeature(
  validSignalStats,
  "distance_m",
  ["ndssi_mean"]
).setChartType("LineChart").setOptions({
  title: "February 2023 NDSSI by Distance Downstream",
  hAxis: {title: "Distance from impact pool (m)"},
  vAxis: {title: "Mean NDSSI"},
  lineWidth: 2,
  pointSize: 5,
  colors: ["#1565C0"]
});
print(ndssiDistanceChart);

print("Script ready. Check the polygon table first; if most polygons have zero or only a few valid pixels, the geometry size or month choice is too constrained for this analysis.");
