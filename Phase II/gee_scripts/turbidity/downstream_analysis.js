/*
  FILE: Phase II/gee_scripts/turbidity/downstream_analysis.js
  PURPOSE: Quick feasibility check for valid Sentinel-2 pixels plus uploaded
           Planet Labs scenes, with monthly EGRI / NDSSI / red-band turbidity
           proxy comparisons across user-drawn downstream polygons.

  GEE IMPORTS REQUIRED
  - upstream_control
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

  PLANET LABS ASSET SETUP
  1. Upload each Planet Labs surface reflectance GeoTIFF as an Earth Engine image.
  2. If available, upload each matching UDM / usable-data-mask GeoTIFF as an image.
  3. For split scenes, upload each tile separately and list them in srAssetIds / udmAssetIds.
  4. Replace the asset placeholders in planetScenes below.
  5. Set enabled: true for scenes that have finished uploading.

  Current Rio Santa Rosa downstream re-orthotile package found locally:
  - 2011-02-24 RapidEye-2 tiles: 1645310, 1645311
  - 2012-02-14 RapidEye-4 tiles: 1645310, 1645311
  - 2012-02-26 RapidEye-2 tiles: 1645310, 1645311

  The `RioSantaRosaDownstreamPL_rescene_basic_analytic` package is a separate
  REScene/basic analytic delivery with per-band files. This script is prepared
  for the re-orthotile surface reflectance assets first.

  Band index presets:
  - RapidEye 5-band SR:       [0, 1, 2, 4] = blue, green, red, nir
  - PlanetScope 4-band SR:    [0, 1, 2, 3] = blue, green, red, nir
  - PlanetScope 8-band SR: use the asset band list to set blue/green/red/nir.
*/

// ================= USER SETTINGS =================
var cloudMax = 60;
var sentinelScaleMeters = 10;
var planetScaleMeters = 5;
var analysisScaleMeters = planetScaleMeters;
var upstreamReferenceDistanceMeters = -500;
var analysisMonths = [
  {label: "2023-02", title: "February 2023", start: ee.Date("2023-02-01")},
  {label: "2024-01", title: "January 2024", start: ee.Date("2024-01-01")}
];

var includePlanetLabsTiles = false;
var planetScenes = [
  {
    enabled: false,
    label: "Rio Santa Rosa downstream RapidEye-2 2011-02-24",
    sensor: "RapidEye-2",
    acquisitionDate: "2011-02-24",
    srAssetIds: [
      "projects/metalminingpersonalcopy/assets/1645310_2011-02-24_RE2_3A_Analytic_SR_clip",
      "projects/metalminingpersonalcopy/assets/1645311_2011-02-24_RE2_3A_Analytic_SR_clip"
    ],
    udmAssetIds: [
      "projects/metalminingpersonalcopy/assets/1645310_2011-02-24_RE2_3A_udm_clip",
      "projects/metalminingpersonalcopy/assets/1645311_2011-02-24_RE2_3A_udm_clip"
    ],
    bandIndexes: [0, 1, 2, 4]
  },
  {
    enabled: false,
    label: "Rio Santa Rosa downstream RapidEye-4 2012-02-14",
    sensor: "RapidEye-4",
    acquisitionDate: "2012-02-14",
    srAssetIds: [
      "projects/metalminingpersonalcopy/assets/1645310_2012-02-14_RE4_3A_Analytic_SR_clip",
      "projects/metalminingpersonalcopy/assets/1645311_2012-02-14_RE4_3A_Analytic_SR_clip"
    ],
    udmAssetIds: [
      "projects/metalminingpersonalcopy/assets/1645310_2012-02-14_RE4_3A_udm_clip",
      "projects/metalminingpersonalcopy/assets/1645311_2012-02-14_RE4_3A_udm_clip"
    ],
    bandIndexes: [0, 1, 2, 4]
  },
  {
    enabled: false,
    label: "Rio Santa Rosa downstream RapidEye-2 2012-02-26",
    sensor: "RapidEye-2",
    acquisitionDate: "2012-02-26",
    srAssetIds: [
      "projects/metalminingpersonalcopy/assets/1645310_2012-02-26_RE2_3A_Analytic_SR_clip",
      "projects/metalminingpersonalcopy/assets/1645311_2012-02-26_RE2_3A_Analytic_SR_clip"
    ],
    udmAssetIds: [
      "projects/metalminingpersonalcopy/assets/1645310_2012-02-26_RE2_3A_udm_clip",
      "projects/metalminingpersonalcopy/assets/1645311_2012-02-26_RE2_3A_udm_clip"
    ],
    bandIndexes: [0, 1, 2, 4]
  },
  {
    enabled: false,
    label: "Optional extra Planet scene",
    sensor: "Planet Labs",
    acquisitionDate: "YYYY-MM-DD",
    srAssetIds: [
      "projects/YOUR_PROJECT/assets/Downstream_Planet_SR_scene_4"
    ],
    udmAssetIds: [
      "projects/YOUR_PROJECT/assets/Downstream_Planet_UDM_scene_4"
    ],
    bandIndexes: [0, 1, 2, 3]
  }
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
  return img.updateMask(bad.not())
    .select(["blue", "green", "red", "nir"])
    .set({
      sensor: "Sentinel-2",
      source: "COPERNICUS/S2_SR_HARMONIZED",
      analysis_scale_m: sentinelScaleMeters
    })
    .copyProperties(img, ["system:time_start", "CLOUDY_PIXEL_PERCENTAGE"]);
}

function scaleAndSelectS2(img) {
  var reflectance = img
    .select(["B2", "B3", "B4", "B8"], ["blue", "green", "red", "nir"])
    .multiply(0.0001);

  return reflectance
    .addBands(img.select(["QA60", "SCL"]))
    .copyProperties(img, ["system:time_start", "CLOUDY_PIXEL_PERCENTAGE"]);
}

function addIndices(img) {
  var ndssi = img.expression(
    "(BLUE - NIR) / (BLUE + NIR)",
    {
      BLUE: img.select("blue"),
      NIR: img.select("nir")
    }
  ).rename("NDSSI");

  var egri = img.expression(
    "GREEN / RED",
    {
      GREEN: img.select("green"),
      RED: img.select("red").max(0.0001)
    }
  ).rename("EGRI");

  var redTurbidityProxy = img.select("red").rename("RED_TURBIDITY_PROXY");

  var hossainRedTurbidity = img.expression(
    "2677.2 * pow(RED, 1.8562)",
    {
      RED: img.select("red").max(0)
    }
  ).rename("HOSSAIN_RED_NTU_PROXY");

  return img.addBands([ndssi, egri, redTurbidityProxy, hossainRedTurbidity]);
}

function makeEmptyAnalysisImage() {
  return ee.Image.constant([0, 0, 0, 0, 0, 0, 0, 0])
    .rename([
      "blue",
      "green",
      "red",
      "nir",
      "NDSSI",
      "EGRI",
      "RED_TURBIDITY_PROXY",
      "HOSSAIN_RED_NTU_PROXY"
    ])
    .clip(aoi)
    .selfMask();
}

function toJsArray(value) {
  if (!value) {
    return [];
  }
  return Array.isArray(value) ? value : [value];
}

function makePlanetLabsImage(scene) {
  var sceneDate = ee.Date(scene.acquisitionDate);
  var srAssetIds = toJsArray(scene.srAssetIds || scene.srAssetId);
  var udmAssetIds = toJsArray(scene.udmAssetIds || scene.udmAssetId);

  var srCollection = ee.ImageCollection(srAssetIds.map(function(assetId) {
    return ee.Image(assetId)
      .select(scene.bandIndexes, ["blue", "green", "red", "nir"])
      .multiply(0.0001);
  }));

  var sr = srCollection.mosaic();

  var positiveReflectanceMask = sr.select("blue").gt(0)
    .and(sr.select("green").gt(0))
    .and(sr.select("red").gt(0))
    .and(sr.select("nir").gt(0));

  var validMask = positiveReflectanceMask;
  if (udmAssetIds.length > 0) {
    var udmCollection = ee.ImageCollection(udmAssetIds.map(function(assetId) {
      return ee.Image(assetId).select([0], ["UDM"]);
    }));
    var udm = udmCollection.mosaic();
    validMask = validMask.and(udm.eq(0));
  }

  return addIndices(sr.updateMask(validMask))
    .set({
      "system:time_start": sceneDate.millis(),
      date: sceneDate.format("YYYY-MM-dd"),
      scene_label: scene.label,
      sensor: scene.sensor,
      source: "Planet Labs uploaded asset",
      sr_asset_ids: srAssetIds.join(","),
      udm_asset_ids: udmAssetIds.join(","),
      tile_count: srAssetIds.length,
      analysis_scale_m: planetScaleMeters
    });
}

function summarizePolygon(feature, imageCollection, monthImage, imageCount) {
  feature = ee.Feature(feature);
  var geom = feature.geometry();
  var areaSqm = geom.area(1);

  var validMask = monthImage.select("red").mask().clip(geom);
  var validPx = ee.Number(ee.Algorithms.If(
    imageCount.gt(0),
    validMask.reduceRegion({
      reducer: ee.Reducer.count(),
      geometry: geom,
      scale: analysisScaleMeters,
      bestEffort: true,
      maxPixels: 1e8
    }).get("red"),
    0
  ));

  var perImageStats = ee.FeatureCollection(ee.Algorithms.If(
    imageCount.gt(0),
    imageCollection.map(function(img) {
      var imageValidMask = img.select("red").mask().clip(geom);
      var imageValidPx = ee.Number(imageValidMask.reduceRegion({
        reducer: ee.Reducer.count(),
        geometry: geom,
        scale: analysisScaleMeters,
        bestEffort: true,
        maxPixels: 1e8
      }).get("red"));

      var imageStats = ee.Dictionary(ee.Algorithms.If(
        imageValidPx.gt(0),
        img.select(["NDSSI", "EGRI", "RED_TURBIDITY_PROXY", "HOSSAIN_RED_NTU_PROXY"])
          .updateMask(imageValidMask)
          .reduceRegion({
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
        ndssi_mean: imageStats.get("NDSSI"),
        egri_mean: imageStats.get("EGRI"),
        red_turbidity_proxy_mean: imageStats.get("RED_TURBIDITY_PROXY"),
        hossain_red_ntu_proxy_mean: imageStats.get("HOSSAIN_RED_NTU_PROXY")
      });
    }),
    ee.FeatureCollection([])
  ));

  var validPerImageStats = perImageStats.filter(ee.Filter.gt("valid_px", 0));
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

  var redTurbidityProxyMedian = ee.Algorithms.If(
    validPerImageStats.size().gt(0),
    ee.Array(validPerImageStats.aggregate_array("red_turbidity_proxy_mean"))
      .reduce(ee.Reducer.median(), [0])
      .get([0]),
    null
  );

  var hossainRedNtuProxyMedian = ee.Algorithms.If(
    validPerImageStats.size().gt(0),
    ee.Array(validPerImageStats.aggregate_array("hossain_red_ntu_proxy_mean"))
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
    red_turbidity_proxy_mean: redTurbidityProxyMedian,
    hossain_red_ntu_proxy_mean: hossainRedNtuProxyMedian,
    per_image_valid_count: validPerImageStats.size()
  });
}

function summarizePolygonForSingleImage(feature, image, imageLabel, imageCount) {
  feature = ee.Feature(feature);
  var geom = feature.geometry();
  var validMask = image.select("red").mask().clip(geom);

  var validPx = ee.Number(ee.Algorithms.If(
    imageCount.gt(0),
    validMask.reduceRegion({
      reducer: ee.Reducer.count(),
      geometry: geom,
      scale: analysisScaleMeters,
      bestEffort: true,
      maxPixels: 1e8
    }).get("red"),
    0
  ));

  var stats = ee.Dictionary(ee.Algorithms.If(
    validPx.gt(0),
    image.select(["NDSSI", "EGRI", "RED_TURBIDITY_PROXY", "HOSSAIN_RED_NTU_PROXY"])
      .updateMask(validMask)
      .reduceRegion({
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
    ndssi_mean: stats.get("NDSSI"),
    egri_mean: stats.get("EGRI"),
    red_turbidity_proxy_mean: stats.get("RED_TURBIDITY_PROXY"),
    hossain_red_ntu_proxy_mean: stats.get("HOSSAIN_RED_NTU_PROXY"),
    image_label: imageLabel
  });
}

// ================= INPUTS =================
var polygons = ee.FeatureCollection([
  ee.Feature(toGeometry(upstream_control), {
    polygon_id: "upstream_control",
    distance_m: upstreamReferenceDistanceMeters
  }),
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
var activePlanetScenes = includePlanetLabsTiles ? planetScenes.filter(function(scene) {
  return scene.enabled;
}) : [];
var planetCollection = ee.ImageCollection(activePlanetScenes.map(makePlanetLabsImage));

print("Sentinel-2 scale (m)", sentinelScaleMeters);
print("Planet Labs scale (m)", planetScaleMeters);
print("Analysis summary scale (m)", analysisScaleMeters);
print("Planet Labs enabled", includePlanetLabsTiles);
print("Active Planet Labs scene count", activePlanetScenes.length);
print("Planet Labs scene config", activePlanetScenes);
print("Upstream reference distance (m)", upstreamReferenceDistanceMeters);
print("NDSSI formula", "(Blue - NIR) / (Blue + NIR)");
print("EGRI formula", "Green / Red");
print("Red turbidity proxy", "Red surface reflectance");
print("Hossain red-band proxy", "2677.2 * pow(red, 1.8562); exploratory, not locally calibrated");

function buildMonthSummary(config) {
  var monthStart = config.start;
  var monthEnd = monthStart.advance(1, "month");

  var s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
    .filterBounds(aoi)
    .filterDate(monthStart, monthEnd)
    .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", cloudMax))
    .map(scaleAndSelectS2)
    .map(maskS2)
    .map(addIndices);

  var planetMonth = planetCollection
    .filterBounds(aoi)
    .filterDate(monthStart, monthEnd);
  var analysisImages = s2.merge(planetMonth).sort("system:time_start");

  var imageCount = analysisImages.size();
  var sortedImages = analysisImages.sort("system:time_start");
  var secondImage = ee.Image(ee.Algorithms.If(
    imageCount.gt(1),
    sortedImages.toList(2).get(1),
    ee.Algorithms.If(
      imageCount.gt(0),
      sortedImages.first(),
      makeEmptyAnalysisImage()
    )
  ));
  var secondImageLabel = ee.String(ee.Algorithms.If(
    imageCount.gt(1),
    ee.String(secondImage.get("sensor")).cat("_").cat(ee.Date(secondImage.get("system:time_start")).format("YYYY-MM-dd HH:mm")),
    ee.Algorithms.If(
      imageCount.gt(0),
      ee.String("fallback_first_image_")
        .cat(ee.String(secondImage.get("sensor")))
        .cat("_")
        .cat(ee.Date(secondImage.get("system:time_start")).format("YYYY-MM-dd HH:mm")),
      "no_valid_image"
    )
  ));
  var monthImage = ee.Image(ee.Algorithms.If(
    imageCount.gt(0),
    analysisImages.median().clip(aoi),
    makeEmptyAnalysisImage()
  ));

  var polygonStats = polygons.map(function(feature) {
    return summarizePolygon(feature, analysisImages, monthImage, imageCount);
  }).sort("distance_m");

  var validSignalStats = polygonStats.filter(ee.Filter.gt("valid_px", 0));
  var secondImagePolygonStats = polygons.map(function(feature) {
    return summarizePolygonForSingleImage(feature, secondImage, secondImageLabel, imageCount);
  }).sort("distance_m");
  var secondImageValidSignalStats = secondImagePolygonStats.filter(ee.Filter.gt("valid_px", 0));

  return {
    title: config.title,
    label: config.label,
    monthStart: monthStart,
    monthEnd: monthEnd,
    imageCount: imageCount,
    sentinelImageCount: s2.size(),
    planetImageCount: planetMonth.size(),
    secondImage: secondImage,
    secondImageLabel: secondImageLabel,
    monthImage: monthImage,
    polygonStats: polygonStats,
    validSignalStats: validSignalStats,
    secondImagePolygonStats: secondImagePolygonStats,
    secondImageValidSignalStats: secondImageValidSignalStats
  };
}

// ================= OUTPUTS =================
Map.setOptions("SATELLITE");
Map.centerObject(aoi, 15);
Map.addLayer(polygons, {color: "FFB300"}, "Downstream polygons", true);
Map.addLayer(toGeometry(upstream_control), {color: "26A69A"}, "upstream_control", true);
Map.addLayer(toGeometry(impact_pool), {color: "FFFF00"}, "impact_pool", true);
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
  return summary.secondImageValidSignalStats.map(function(feature) {
    return ee.Feature(feature).set({
      month_label: summary.title,
      month_key: summary.label,
      image_label: summary.secondImageLabel
    });
  });
})).flatten();
var planetSceneSignalStats = ee.FeatureCollection(activePlanetScenes.map(function(scene) {
  var sceneImage = makePlanetLabsImage(scene).clip(aoi);
  var sceneLabel = scene.label;

  return polygons.map(function(feature) {
    return summarizePolygonForSingleImage(feature, sceneImage, sceneLabel, ee.Number(1))
      .set({
        scene_label: sceneLabel,
        acquisition_date: scene.acquisitionDate,
        sensor_label: scene.sensor
      });
  }).filter(ee.Filter.gt("valid_px", 0));
})).flatten();

monthSummaries.forEach(function(summary, idx) {
  Map.addLayer(
    summary.monthImage.select(["red", "green", "blue"]),
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

  Map.addLayer(
    summary.monthImage.select("RED_TURBIDITY_PROXY"),
    {min: 0.01, max: 0.18, palette: ["#F7FCF0", "#FEE08B", "#D73027"]},
    summary.title + " red turbidity proxy",
    false
  );

  print(summary.title + " window", summary.monthStart.format("YYYY-MM-dd"), summary.monthEnd.format("YYYY-MM-dd"));
  print(summary.title + " Sentinel-2 image count after filtering", summary.sentinelImageCount);
  print(summary.title + " Planet Labs image count", summary.planetImageCount);
  print(summary.title + " combined image count", summary.imageCount);
  print(summary.title + " second image used for single-image charts", summary.secondImageLabel);
  print(summary.title + " valid pixel feasibility by polygon", summary.polygonStats);
  print(summary.title + " polygons with valid pixels", summary.polygonStats.filter(ee.Filter.gt("valid_px", 0)));
  print(summary.title + " polygons with zero valid pixels", summary.polygonStats.filter(ee.Filter.eq("valid_px", 0)));
  print(summary.title + " polygon signal summary for comparison", summary.validSignalStats);
  print(summary.title + " second-image polygon signal summary", summary.secondImageValidSignalStats);

  print(ui.Chart.feature.byFeature(
    summary.polygonStats,
    "polygon_id",
    ["valid_px"]
  ).setChartType("ColumnChart").setOptions({
    title: summary.title + " Valid Pixels by Polygon",
    hAxis: {title: "Polygon"},
    vAxis: {title: "Valid pixels at analysis scale"},
    colors: ["#FB8C00"]
  }));
});

print(ui.Chart.feature.groups(
  combinedValidSignalStats,
  "distance_m",
  "egri_mean",
  "month_label"
).setChartType("LineChart").setOptions({
  title: "EGRI by Distance Downstream",
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
  title: "NDSSI by Distance Downstream",
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
  title: "EGRI by Distance Downstream (Second Image of Month)",
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
  title: "NDSSI by Distance Downstream (Second Image of Month)",
  hAxis: {title: "Distance from impact pool (m)"},
  vAxis: {title: "Mean NDSSI"},
  lineWidth: 2,
  pointSize: 5,
  colors: ["#0D47A1", "#4A148C"]
}));

print(ui.Chart.feature.groups(
  combinedValidSignalStats,
  "distance_m",
  "red_turbidity_proxy_mean",
  "month_label"
).setChartType("LineChart").setOptions({
  title: "Red Turbidity Proxy by Distance Downstream",
  hAxis: {title: "Distance from impact pool (m)"},
  vAxis: {title: "Mean red reflectance"},
  lineWidth: 2,
  pointSize: 5,
  colors: ["#EF6C00", "#AD1457"]
}));

print(ui.Chart.feature.groups(
  combinedValidSignalStats,
  "distance_m",
  "hossain_red_ntu_proxy_mean",
  "month_label"
).setChartType("LineChart").setOptions({
  title: "Hossain Red-Band NTU Proxy by Distance Downstream",
  hAxis: {title: "Distance from impact pool (m)"},
  vAxis: {title: "Exploratory NTU proxy"},
  lineWidth: 2,
  pointSize: 5,
  colors: ["#E65100", "#880E4F"]
}));

print("Planet Labs per-scene downstream signal summary", planetSceneSignalStats);

print(ui.Chart.feature.groups(
  planetSceneSignalStats,
  "distance_m",
  "egri_mean",
  "scene_label"
).setChartType("LineChart").setOptions({
  title: "Planet Labs EGRI by Distance Downstream",
  hAxis: {title: "Distance from impact pool (m)"},
  vAxis: {title: "Mean EGRI"},
  lineWidth: 2,
  pointSize: 5,
  colors: ["#1B5E20", "#2E7D32", "#558B2F", "#33691E"]
}));

print(ui.Chart.feature.groups(
  planetSceneSignalStats,
  "distance_m",
  "ndssi_mean",
  "scene_label"
).setChartType("LineChart").setOptions({
  title: "Planet Labs NDSSI by Distance Downstream",
  hAxis: {title: "Distance from impact pool (m)"},
  vAxis: {title: "Mean NDSSI"},
  lineWidth: 2,
  pointSize: 5,
  colors: ["#0D47A1", "#1565C0", "#3949AB", "#6A1B9A"]
}));

print(ui.Chart.feature.groups(
  planetSceneSignalStats,
  "distance_m",
  "hossain_red_ntu_proxy_mean",
  "scene_label"
).setChartType("LineChart").setOptions({
  title: "Planet Labs Hossain Red-Band Proxy by Distance Downstream",
  hAxis: {title: "Distance from impact pool (m)"},
  vAxis: {title: "Exploratory NTU proxy"},
  lineWidth: 2,
  pointSize: 5,
  colors: ["#BF360C", "#E65100", "#F57C00", "#AD1457"]
}));

activePlanetScenes.forEach(function(scene, index) {
  var sceneImage = makePlanetLabsImage(scene).clip(aoi);
  Map.addLayer(
    sceneImage.select(["red", "green", "blue"]),
    {min: 0.02, max: 0.35, gamma: 1.2},
    "Planet Labs RGB " + scene.label,
    index === activePlanetScenes.length - 1
  );
  Map.addLayer(
    sceneImage.select("EGRI"),
    {min: 0.5, max: 2.0, palette: ["#7F0000", "#FDD49E", "#238B45"]},
    "Planet Labs EGRI " + scene.label,
    false
  );
  Map.addLayer(
    sceneImage.select("NDSSI"),
    {min: -0.8, max: 0.4, palette: ["#54278F", "#2B8CBE", "#F7FCF0"]},
    "Planet Labs NDSSI " + scene.label,
    false
  );
});

print("Script ready. Fill the four Planet Labs scene placeholders, set enabled: true for uploaded scenes, and set includePlanetLabsTiles: true to include them in downstream polygon summaries.");
