/*
  FILE: Phase I/gee_scripts/turbidity/san_sebastian_downstream_polygon_valid_pixels_feb2023.js
  PURPOSE: Quick feasibility check for valid Sentinel-2 pixels plus monthly
           EGRI / NDSSI distance comparisons across user-drawn downstream
           polygons for February 2023 and January 2024.

  METHODOLOGY NOTE
  - Monthly summaries and single-image comparisons only use images where
    all but one polygon have valid pixels, at minimum. Partial-coverage
    scenes below that threshold are excluded.

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
var allowedMissingPolygons = 1;
var analysisMonths = [
  {label: "2023-02", title: "February 2023", start: ee.Date("2023-02-01")},
  {label: "2024-01", title: "January 2024", start: ee.Date("2024-01-01")}
];

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

function dictValueOrNull(dict, key) {
  dict = ee.Dictionary(dict);
  return ee.Algorithms.If(dict.contains(key), dict.get(key), null);
}

function countMaskedPixels(image, geom) {
  var count = image.reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: geom,
    scale: analysisScaleMeters,
    bestEffort: true,
    maxPixels: 1e8
  }).get("B4");

  return ee.Number(ee.Algorithms.If(ee.Algorithms.IsEqual(count, null), 0, count));
}

function summarizePolygon(feature, imageCollection, monthImage, imageCount) {
  feature = ee.Feature(feature);
  var geom = feature.geometry();
  var areaSqm = geom.area(1);

  var validMask = monthImage.select("B4").mask().clip(geom);
  var validPx = ee.Number(ee.Algorithms.If(
    imageCount.gt(0),
    countMaskedPixels(validMask, geom),
    0
  ));

  var perImageStats = ee.FeatureCollection(ee.Algorithms.If(
    imageCount.gt(0),
    imageCollection.map(function(img) {
      var imageValidMask = img.select("B4").mask().clip(geom);
      var imageValidPx = countMaskedPixels(imageValidMask, geom);

      var imageStats = ee.Dictionary(ee.Algorithms.If(
        imageValidPx.gt(0),
        img.select(["NDSSI", "EGRI"]).updateMask(imageValidMask).reduceRegion({
          reducer: ee.Reducer.mean(),
          geometry: geom,
          scale: analysisScaleMeters,
          bestEffort: true,
          maxPixels: 1e8
        }),
        ee.Dictionary({})
      ));

      return ee.Feature(null, {
        valid_px: imageValidPx,
        ndssi_mean: dictValueOrNull(imageStats, "NDSSI"),
        egri_mean: dictValueOrNull(imageStats, "EGRI")
      });
    }),
    ee.FeatureCollection([])
  ));

  var validPerImageStats = perImageStats
    .filter(ee.Filter.gt("valid_px", 0))
    .filter(ee.Filter.notNull(["ndssi_mean", "egri_mean"]));
  var ndssiMedian = ee.Algorithms.If(
    validPerImageStats.size().gt(0),
    ee.Array(validPerImageStats.aggregate_array("ndssi_mean"))
      .reduce(ee.Reducer.median(), [0])
      .get([0]),
    null
  );

  var egriMedian = ee.Algorithms.If(
    validPerImageStats.size().gt(0),
    ee.Array(validPerImageStats.aggregate_array("egri_mean"))
      .reduce(ee.Reducer.median(), [0])
      .get([0]),
    null
  );

  return feature.set({
    polygon_area_sqm: areaSqm,
    image_count: imageCount,
    valid_px: validPx,
    has_valid_pixels: validPx.gt(0),
    ndssi_mean: ndssiMedian,
    egri_mean: egriMedian,
    per_image_valid_count: validPerImageStats.size()
  });
}

function summarizePolygonForSingleImage(feature, image, imageLabel, imageCount) {
  feature = ee.Feature(feature);
  var geom = feature.geometry();
  var validMask = image.select("B4").mask().clip(geom);

  var validPx = ee.Number(ee.Algorithms.If(
    imageCount.gt(0),
    countMaskedPixels(validMask, geom),
    0
  ));

  var stats = ee.Dictionary(ee.Algorithms.If(
    validPx.gt(0),
    image.select(["NDSSI", "EGRI"]).updateMask(validMask).reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: geom,
      scale: analysisScaleMeters,
      bestEffort: true,
      maxPixels: 1e8
    }),
    ee.Dictionary({})
  ));

  return feature.set({
    valid_px: validPx,
    has_valid_pixels: validPx.gt(0),
    ndssi_mean: dictValueOrNull(stats, "NDSSI"),
    egri_mean: dictValueOrNull(stats, "EGRI"),
    image_label: imageLabel
  });
}

function imageHasAllPolygonsValid(image) {
  var validityChecks = polygons.map(function(feature) {
    var geom = ee.Feature(feature).geometry();
    var validPx = countMaskedPixels(image.select("B4").mask().clip(geom), geom);

    return ee.Feature(null, {has_valid: validPx.gt(0)});
  });

  var validCount = ee.FeatureCollection(validityChecks)
    .filter(ee.Filter.eq("has_valid", true))
    .size();

  var minValidPolygons = polygons.size().subtract(allowedMissingPolygons);
  return image.set({
    all_polygons_valid: validCount.eq(polygons.size()),
    enough_polygons_valid: validCount.gte(minValidPolygons),
    valid_polygon_count: validCount,
    missing_polygon_count: polygons.size().subtract(validCount)
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

function buildMonthSummary(config) {
  var monthStart = config.start;
  var monthEnd = monthStart.advance(1, "month");

  var s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
    .filterBounds(aoi)
    .filterDate(monthStart, monthEnd)
    .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", cloudMax))
    .map(scaleAndSelectS2)
    .map(maskS2)
    .map(addIndices)
    .map(imageHasAllPolygonsValid);

  var imageCount = s2.size();
  var coverageDiagnostics = ee.FeatureCollection(s2.map(function(img) {
    return ee.Feature(null, {
      system_time_start: img.get("system:time_start"),
      image_time: ee.Date(img.get("system:time_start")).format("YYYY-MM-dd HH:mm"),
      valid_polygon_count: img.get("valid_polygon_count"),
      missing_polygon_count: img.get("missing_polygon_count"),
      enough_polygons_valid: img.get("enough_polygons_valid")
    });
  })).sort("system_time_start");
  var fullyValidImages = s2.filter(ee.Filter.eq("enough_polygons_valid", true));
  var fullyValidImageCount = fullyValidImages.size();
  var firstImage = ee.Image(ee.Algorithms.If(
    fullyValidImageCount.gt(0),
    fullyValidImages.sort("system:time_start").first(),
    ee.Image.constant([0, 0, 0, 0, 0, 0]).rename(["B2", "B3", "B4", "B8", "NDSSI", "EGRI"]).clip(aoi).selfMask()
  ));
  var firstImageLabel = ee.String(ee.Algorithms.If(
    fullyValidImageCount.gt(0),
    ee.Date(firstImage.get("system:time_start")).format("YYYY-MM-dd HH:mm"),
    "no_image_meeting_polygon_threshold"
  ));
  var monthImage = ee.Image(ee.Algorithms.If(
    fullyValidImageCount.gt(0),
    fullyValidImages.median().clip(aoi),
    ee.Image.constant([0, 0, 0, 0, 0, 0]).rename(["B2", "B3", "B4", "B8", "NDSSI", "EGRI"]).clip(aoi).selfMask()
  ));

  var polygonStats = polygons.map(function(feature) {
    return summarizePolygon(feature, fullyValidImages, monthImage, fullyValidImageCount);
  }).sort("distance_m");

  var validSignalStats = polygonStats.filter(ee.Filter.gt("valid_px", 0));
  var firstImagePolygonStats = polygons.map(function(feature) {
    return summarizePolygonForSingleImage(feature, firstImage, firstImageLabel, fullyValidImageCount);
  }).sort("distance_m");
  var firstImageValidSignalStats = firstImagePolygonStats.filter(ee.Filter.gt("valid_px", 0));

  return {
    title: config.title,
    label: config.label,
    monthStart: monthStart,
    monthEnd: monthEnd,
    imageCount: imageCount,
    coverageDiagnostics: coverageDiagnostics,
    fullyValidImageCount: fullyValidImageCount,
    firstImage: firstImage,
    firstImageLabel: firstImageLabel,
    monthImage: monthImage,
    polygonStats: polygonStats,
    validSignalStats: validSignalStats,
    firstImagePolygonStats: firstImagePolygonStats,
    firstImageValidSignalStats: firstImageValidSignalStats
  };
}

// ================= OUTPUTS =================
Map.setOptions("SATELLITE");
Map.centerObject(aoi, 15);
Map.addLayer(polygons, {color: "FFB300"}, "Downstream polygons", true);
var monthSummaries = analysisMonths.map(buildMonthSummary);
var combinedValidSignalStats = ee.FeatureCollection(monthSummaries.map(function(summary) {
  return summary.validSignalStats.map(function(feature) {
    return ee.Feature(feature).set({
      month_label: summary.title,
      month_key: summary.label
    });
  });
})).flatten();
var combinedFirstImageSignalStats = ee.FeatureCollection(monthSummaries.map(function(summary) {
  return summary.firstImageValidSignalStats.map(function(feature) {
    return ee.Feature(feature).set({
      month_label: summary.title,
      month_key: summary.label,
      image_label: summary.firstImageLabel
    });
  });
})).flatten();

monthSummaries.forEach(function(summary, idx) {
  Map.addLayer(
    summary.monthImage.select(["B4", "B3", "B2"]),
    {min: 0.02, max: 0.3},
    summary.title + " RGB",
    idx === 0
  );

  Map.addLayer(
    summary.monthImage.select("NDSSI"),
    {min: -0.4, max: 0.4, palette: ["#7B3294", "#F7F7F7", "#008837"]},
    summary.title + " NDSSI",
    false
  );

  print(summary.title + " window", summary.monthStart.format("YYYY-MM-dd"), summary.monthEnd.format("YYYY-MM-dd"));
  print(summary.title + " Sentinel-2 image count after filtering", summary.imageCount);
  print("Allowed missing polygons", allowedMissingPolygons);
  print(summary.title + " image coverage diagnostics", summary.coverageDiagnostics);
  print(summary.title + " images meeting polygon coverage threshold", summary.fullyValidImageCount);
  print(summary.title + " first image meeting polygon coverage threshold used for single-image charts", summary.firstImageLabel);
  print(summary.title + " valid pixel feasibility by polygon (coverage-threshold subset only)", summary.polygonStats);
  print(summary.title + " polygons with valid pixels (coverage-threshold subset only)", summary.polygonStats.filter(ee.Filter.gt("valid_px", 0)));
  print(summary.title + " polygons with zero valid pixels (coverage-threshold subset only)", summary.polygonStats.filter(ee.Filter.eq("valid_px", 0)));
  print(summary.title + " polygon signal summary for comparison (coverage-threshold subset only)", summary.validSignalStats);
  print(summary.title + " first-image polygon signal summary", summary.firstImageValidSignalStats);

  print(ui.Chart.feature.byFeature(
    summary.polygonStats,
    "polygon_id",
    ["valid_px"]
  ).setChartType("ColumnChart").setOptions({
    title: summary.title + " Valid Pixels by Polygon (Coverage Threshold Applied)",
    hAxis: {title: "Polygon"},
    vAxis: {title: "Valid pixels at 10 m"},
    colors: ["#FB8C00"]
  }));
});

print(ui.Chart.feature.groups(
  combinedValidSignalStats,
  "distance_m",
  "egri_mean",
  "month_label"
).setChartType("LineChart").setOptions({
  title: "EGRI by Distance Downstream (Coverage Threshold Applied)",
  hAxis: {title: "Distance from impact pool (m)"},
  vAxis: {title: "Mean EGRI"},
  lineWidth: 2,
  pointSize: 5,
  colors: ["#2E7D32", "#C62828"]
}));

print(ui.Chart.feature.groups(
  combinedValidSignalStats,
  "distance_m",
  "ndssi_mean",
  "month_label"
).setChartType("LineChart").setOptions({
  title: "NDSSI by Distance Downstream (Coverage Threshold Applied)",
  hAxis: {title: "Distance from impact pool (m)"},
  vAxis: {title: "Mean NDSSI"},
  lineWidth: 2,
  pointSize: 5,
  colors: ["#1565C0", "#6A1B9A"]
}));

print(ui.Chart.feature.groups(
  combinedFirstImageSignalStats,
  "distance_m",
  "egri_mean",
  "image_label"
).setChartType("LineChart").setOptions({
  title: "EGRI by Distance Downstream (First Valid Image of Month)",
  hAxis: {title: "Distance from impact pool (m)"},
  vAxis: {title: "Mean EGRI"},
  lineWidth: 2,
  pointSize: 5,
  colors: ["#1B5E20", "#B71C1C"]
}));

print(ui.Chart.feature.groups(
  combinedFirstImageSignalStats,
  "distance_m",
  "ndssi_mean",
  "image_label"
).setChartType("LineChart").setOptions({
  title: "NDSSI by Distance Downstream (First Valid Image of Month)",
  hAxis: {title: "Distance from impact pool (m)"},
  vAxis: {title: "Mean NDSSI"},
  lineWidth: 2,
  pointSize: 5,
  colors: ["#0D47A1", "#4A148C"]
}));

print("Script ready. Monthly summaries and first-image charts now require all but one polygon to have valid coverage. If you meant February 2024 instead of January 2024, change analysisMonths accordingly.");
