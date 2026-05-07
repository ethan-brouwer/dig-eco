/*
  FILE: Phase II/gee_scripts/turbidity/downstream_planet_s2_unified_export.js
  PURPOSE: Export aligned downstream polygon summary CSVs for Sentinel-2 and
           uploaded Planet Labs / RapidEye scenes, including NDTI and raw band
           means, with a shared schema designed for downstream comparison in
           Python or Jupyter.

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

  PLANET / RAPIDEYE ASSET REQUIREMENTS
  - Upload each surface reflectance tile as an ee.Image.
  - Upload each matching UDM tile as an ee.Image.
  - Replace the placeholder asset IDs in planetScenes below.
  - The preferred assets are the re-clipped scenes that include the full
    downstream polygon stack plus the upstream control.

  DESIGN NOTES
  - This script is export-oriented only. It intentionally omits charting and
    map visualization so the logic stays clean and the outputs stay aligned.
  - Both exports use the same column names where possible.
  - Raw band means are included so additional indices, including NDTI-derived
    checks, can be recomputed locally if needed.
*/

// ================= USER SETTINGS =================
var exportFolder = "EarthEngine";
var exportSentinelCsv = true;
var exportPlanetCsv = true;
var exportCombinedCsv = false;

var analysisScaleSentinelMeters = 10;
var analysisScalePlanetMeters = 5;
var aoiBufferMeters = 250;
var upstreamReferenceDistanceMeters = -500;
var minValidPixels = 3;
var maxSentinelScenesPerMonth = 3;
var summaryTileScale = 4;
var reduceRegionMaxPixels = 1e9;
var requirePlanetUdmClearMask = true;

var sentinelStartYear = 2017;
var sentinelEndYear = 2026;
var sentinelMonths = [2, 3]; // February and March
var minSceneClearFractionOverAoi = 0.1;

var includePlanetScenes = true;
var planetScenes = [
  {
    label: "RapidEye-2 2011-02-24",
    sensor: "RapidEye-2",
    acquisitionDate: "2011-02-24",
    srAssetIds: [
      "projects/metalminingpersonalcopy/assets/1645310_2011-02-24_RE2_3A_Analytic_SR_clip_all",
      "projects/metalminingpersonalcopy/assets/1645311_2011-02-24_RE2_3A_Analytic_SR_clip_all"
    ],
    udmAssetIds: [
      "projects/metalminingpersonalcopy/assets/1645310_2011-02-24_RE2_3A_udm_clip_all",
      "projects/metalminingpersonalcopy/assets/1645311_2011-02-24_RE2_3A_udm_clip_all"
    ]
  },
  {
    label: "RapidEye-5 2011-02-27",
    sensor: "RapidEye-5",
    acquisitionDate: "2011-02-27",
    srAssetIds: [
      "projects/metalminingpersonalcopy/assets/1645310_2011-02-27_RE5_3A_Analytic_SR_clip_all",
      "projects/metalminingpersonalcopy/assets/1645311_2011-02-27_RE5_3A_Analytic_SR_clip_all"
    ],
    udmAssetIds: [
      "projects/metalminingpersonalcopy/assets/1645310_2011-02-27_RE5_3A_udm_clip_all",
      "projects/metalminingpersonalcopy/assets/1645311_2011-02-27_RE5_3A_udm_clip_all"
    ]
  },
  {
    label: "RapidEye-4 2012-02-14",
    sensor: "RapidEye-4",
    acquisitionDate: "2012-02-14",
    srAssetIds: [
      "projects/metalminingpersonalcopy/assets/1645310_2012-02-14_RE4_3A_Analytic_SR_clip_all",
      "projects/metalminingpersonalcopy/assets/1645311_2012-02-14_RE4_3A_Analytic_SR_clip_all"
    ],
    udmAssetIds: [
      "projects/metalminingpersonalcopy/assets/1645310_2012-02-14_RE4_3A_udm_clip_all",
      "projects/metalminingpersonalcopy/assets/1645311_2012-02-14_RE4_3A_udm_clip_all"
    ]
  },
  {
    label: "RapidEye-4 2012-03-18",
    sensor: "RapidEye-4",
    acquisitionDate: "2012-03-18",
    srAssetIds: [
      "projects/metalminingpersonalcopy/assets/1645310_2012-03-18_RE4_3A_Analytic_SR_clip_all",
      "projects/metalminingpersonalcopy/assets/1645311_2012-03-18_RE4_3A_Analytic_SR_clip_all"
    ],
    udmAssetIds: [
      "projects/metalminingpersonalcopy/assets/1645310_2012-03-18_RE4_3A_udm_clip_all",
      "projects/metalminingpersonalcopy/assets/1645311_2012-03-18_RE4_3A_udm_clip_all"
    ]
  }
];

var noDataValue = -9999;

// ================= FIELD CONFIGURATION =================
var exportFieldOrder = [
  "export_group",
  "dataset_label",
  "sensor_label",
  "source_label",
  "scene_label",
  "image_label",
  "analysis_date",
  "acquisition_date",
  "analysis_year",
  "analysis_month",
  "month_key",
  "month_label",
  "window_start",
  "window_end_exclusive",
  "image_count",
  "per_image_valid_count",
  "analysis_scale_m",
  "polygon_id",
  "polygon_role",
  "distance_m",
  "polygon_area_sqm",
  "valid_px",
  "has_valid_pixels",
  "qa_flag",
  "blue_mean",
  "green_mean",
  "red_mean",
  "nir_mean",
  "ndssi_mean",
  "egri_mean",
  "ndti_mean",
  "red_turbidity_proxy_mean",
  "hossain_red_ntu_proxy_mean",
  "ndssi_rel_impact_pool",
  "egri_rel_impact_pool",
  "ndti_rel_impact_pool",
  "red_rel_impact_pool",
  "hossain_rel_impact_pool"
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

function safeNumber(x) {
  return ee.Number(ee.Algorithms.If(ee.Algorithms.IsEqual(x, null), noDataValue, x));
}

function safeText(x) {
  return ee.String(ee.Algorithms.If(ee.Algorithms.IsEqual(x, null), "", x));
}

function buildPolygons() {
  var polygonFeatures = [
    ee.Feature(toGeometry(upstream_control), {
      polygon_id: "upstream_control",
      polygon_role: "upstream_control",
      distance_m: upstreamReferenceDistanceMeters
    }),
    ee.Feature(toGeometry(impact_pool), {
      polygon_id: "impact_pool",
      polygon_role: "impact_pool",
      distance_m: 0
    }),
    ee.Feature(toGeometry(Poly500m), {polygon_id: "Poly500m", polygon_role: "downstream", distance_m: 500}),
    ee.Feature(toGeometry(Poly1000m), {polygon_id: "Poly1000m", polygon_role: "downstream", distance_m: 1000}),
    ee.Feature(toGeometry(Poly1500m), {polygon_id: "Poly1500m", polygon_role: "downstream", distance_m: 1500}),
    ee.Feature(toGeometry(Poly2000m), {polygon_id: "Poly2000m", polygon_role: "downstream", distance_m: 2000}),
    ee.Feature(toGeometry(Poly2500m), {polygon_id: "Poly2500m", polygon_role: "downstream", distance_m: 2500}),
    ee.Feature(toGeometry(Poly3000m), {polygon_id: "Poly3000m", polygon_role: "downstream", distance_m: 3000}),
    ee.Feature(toGeometry(Poly3500m), {polygon_id: "Poly3500m", polygon_role: "downstream", distance_m: 3500}),
    ee.Feature(toGeometry(Poly4000m), {polygon_id: "Poly4000m", polygon_role: "downstream", distance_m: 4000}),
    ee.Feature(toGeometry(Poly4500m), {polygon_id: "Poly4500m", polygon_role: "downstream", distance_m: 4500}),
    ee.Feature(toGeometry(Poly5000m), {polygon_id: "Poly5000m", polygon_role: "downstream", distance_m: 5000}),
    ee.Feature(toGeometry(Poly6000m), {polygon_id: "Poly6000m", polygon_role: "downstream", distance_m: 6000}),
    ee.Feature(toGeometry(Poly7000m), {polygon_id: "Poly7000m", polygon_role: "downstream", distance_m: 7000}),
    ee.Feature(toGeometry(Poly8000m), {polygon_id: "Poly8000m", polygon_role: "downstream", distance_m: 8000}),
    ee.Feature(toGeometry(Poly9000m), {polygon_id: "Poly9000m", polygon_role: "downstream", distance_m: 9000}),
    ee.Feature(toGeometry(Poly10000m), {polygon_id: "Poly10000m", polygon_role: "downstream", distance_m: 10000})
  ];
  return ee.FeatureCollection(polygonFeatures).sort("distance_m");
}

function addCoreIndices(img) {
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

  var ndti = img.expression(
    "(RED - GREEN) / (RED + GREEN)",
    {
      RED: img.select("red"),
      GREEN: img.select("green")
    }
  ).rename("NDTI");

  var redTurbidityProxy = img.select("red").rename("RED_TURBIDITY_PROXY");

  var hossainRedTurbidity = img.expression(
    "2677.2 * pow(RED, 1.8562)",
    {
      RED: img.select("red").max(0)
    }
  ).rename("HOSSAIN_RED_NTU_PROXY");

  return img.addBands([
    ndssi,
    egri,
    ndti,
    redTurbidityProxy,
    hossainRedTurbidity
  ]);
}

function attachMatchedProperty(baseFc, lookupFc, groupProperty, lookupProperty, outputProperty) {
  return baseFc.map(function(feature) {
    feature = ee.Feature(feature);
    var match = ee.Feature(
      lookupFc
        .filter(ee.Filter.eq(groupProperty, feature.get(groupProperty)))
        .filter(ee.Filter.eq("polygon_id", feature.get("polygon_id")))
        .first()
    );
    var matchedValue = ee.Algorithms.If(
      ee.Algorithms.IsEqual(match, null),
      null,
      match.get(lookupProperty)
    );
    return feature.set(outputProperty, matchedValue);
  });
}

function normalizeByImpactPool(fc, groupProperty, valueProperty, outputProperty) {
  var groupValues = ee.List(fc.aggregate_array(groupProperty)).distinct();
  return ee.FeatureCollection(groupValues.map(function(groupValue) {
    var groupFc = fc.filter(ee.Filter.eq(groupProperty, groupValue));
    var baselineFeature = ee.Feature(groupFc.filter(ee.Filter.eq("polygon_id", "impact_pool")).first());
    var baselineRaw = baselineFeature.get(valueProperty);
    var baselineValue = ee.Number(baselineRaw);
    var baselineIsValid = ee.Algorithms.If(
      ee.Algorithms.IsEqual(baselineRaw, null),
      false,
      baselineValue.neq(noDataValue).and(baselineValue.neq(0))
    );

    return groupFc.map(function(feature) {
      feature = ee.Feature(feature);
      var valueRaw = feature.get(valueProperty);
      var value = ee.Number(valueRaw);
      var valueIsValid = ee.Algorithms.If(
        ee.Algorithms.IsEqual(valueRaw, null),
        false,
        value.neq(noDataValue)
      );
      var normalizedValue = ee.Algorithms.If(
        baselineIsValid.and(valueIsValid),
        value.divide(baselineValue),
        noDataValue
      );
      return feature.set(outputProperty, normalizedValue);
    });
  })).flatten();
}

function addRelativeFields(fc, groupProperty) {
  var egriNorm = normalizeByImpactPool(fc, groupProperty, "egri_mean", "egri_rel_impact_pool");
  var ndssiNorm = normalizeByImpactPool(fc, groupProperty, "ndssi_mean", "ndssi_rel_impact_pool");
  var ndtiNorm = normalizeByImpactPool(fc, groupProperty, "ndti_mean", "ndti_rel_impact_pool");
  var redNorm = normalizeByImpactPool(fc, groupProperty, "red_turbidity_proxy_mean", "red_rel_impact_pool");
  var hossainNorm = normalizeByImpactPool(fc, groupProperty, "hossain_red_ntu_proxy_mean", "hossain_rel_impact_pool");

  return attachMatchedProperty(
    attachMatchedProperty(
      attachMatchedProperty(
        attachMatchedProperty(
          attachMatchedProperty(
            fc,
            egriNorm,
            groupProperty,
            "egri_rel_impact_pool",
            "egri_rel_impact_pool"
          ),
          ndssiNorm,
          groupProperty,
          "ndssi_rel_impact_pool",
          "ndssi_rel_impact_pool"
        ),
        ndtiNorm,
        groupProperty,
        "ndti_rel_impact_pool",
        "ndti_rel_impact_pool"
      ),
      redNorm,
      groupProperty,
      "red_rel_impact_pool",
      "red_rel_impact_pool"
    ),
    hossainNorm,
    groupProperty,
    "hossain_rel_impact_pool",
    "hossain_rel_impact_pool"
  );
}

function selectExportFields(fc) {
  return fc.map(function(feature) {
    return ee.Feature(feature).select(exportFieldOrder, exportFieldOrder, false);
  });
}

function hasPlaceholderAsset(scene) {
  var srAssetIds = scene.srAssetIds || [];
  var udmAssetIds = scene.udmAssetIds || [];
  var hasBadSr = srAssetIds.some(function(assetId) {
    return String(assetId).indexOf("REPLACE_WITH_") === 0;
  });
  var hasBadUdm = udmAssetIds.some(function(assetId) {
    return String(assetId).indexOf("REPLACE_WITH_") === 0;
  });
  return srAssetIds.length === 0 ||
    udmAssetIds.length === 0 ||
    srAssetIds.length !== udmAssetIds.length ||
    hasBadSr ||
    hasBadUdm ||
    String(scene.acquisitionDate).indexOf("YYYY-") === 0;
}

function filterConfiguredPlanetScenes(sceneList) {
  return sceneList.filter(function(scene) {
    return !hasPlaceholderAsset(scene);
  });
}

function emptyAnalysisImage() {
  return ee.Image.constant([0, 0, 0, 0, 0, 0, 0, 0, 0])
    .rename([
      "blue",
      "green",
      "red",
      "nir",
      "NDSSI",
      "EGRI",
      "NDTI",
      "RED_TURBIDITY_PROXY",
      "HOSSAIN_RED_NTU_PROXY"
    ])
    .selfMask();
}

function summarizeImageToPolygons(image, polygons, config) {
  var metricBandNames = [
    "blue",
    "green",
    "red",
    "nir",
    "NDSSI",
    "EGRI",
    "NDTI",
    "RED_TURBIDITY_PROXY",
    "HOSSAIN_RED_NTU_PROXY"
  ];

  var meanReducer = ee.Reducer.mean();
  var meanStats = ee.FeatureCollection(ee.Algorithms.If(
    ee.Number(config.image_count).gt(0),
    image.select(metricBandNames).reduceRegions({
      collection: polygons,
      reducer: meanReducer,
      scale: config.analysis_scale_m,
      tileScale: summaryTileScale
    }),
    polygons
  ));

  var validPxImage = ee.Image.constant(1)
    .updateMask(image.select("red").mask())
    .rename("valid_px");

  var validPxStats = ee.FeatureCollection(ee.Algorithms.If(
    ee.Number(config.image_count).gt(0),
    validPxImage.reduceRegions({
      collection: polygons,
      reducer: ee.Reducer.count(),
      scale: config.analysis_scale_m,
      tileScale: summaryTileScale
    }),
    polygons.map(function(feature) {
      return ee.Feature(feature).set("count", 0);
    })
  ));

  return attachMatchedProperty(
    meanStats,
    validPxStats,
    "polygon_id",
    "count",
    "valid_px"
  ).map(function(feature) {
    feature = ee.Feature(feature);
    var validPx = ee.Number(ee.Algorithms.If(
      ee.Algorithms.IsEqual(feature.get("valid_px"), null),
      0,
      feature.get("valid_px")
    ));

    return feature.set({
      export_group: config.export_group,
      dataset_label: config.dataset_label,
      sensor_label: config.sensor_label,
      source_label: config.source_label,
      scene_label: config.scene_label,
      image_label: config.image_label,
      analysis_date: config.analysis_date,
      acquisition_date: safeText(config.acquisition_date),
      analysis_year: config.analysis_year,
      analysis_month: config.analysis_month,
      month_key: config.month_key,
      month_label: config.month_label,
      window_start: config.window_start,
      window_end_exclusive: config.window_end_exclusive,
      image_count: config.image_count,
      per_image_valid_count: config.per_image_valid_count,
      analysis_scale_m: config.analysis_scale_m,
      polygon_area_sqm: feature.geometry().area(1),
      valid_px: validPx,
      has_valid_pixels: validPx.gte(minValidPixels),
      qa_flag: ee.Algorithms.If(
        ee.Number(config.image_count).eq(0),
        "no_images",
        ee.Algorithms.If(validPx.gte(minValidPixels), "ok", "low_valid_px")
      ),
      blue_mean: safeNumber(feature.get("blue")),
      green_mean: safeNumber(feature.get("green")),
      red_mean: safeNumber(feature.get("red")),
      nir_mean: safeNumber(feature.get("nir")),
      ndssi_mean: safeNumber(feature.get("NDSSI")),
      egri_mean: safeNumber(feature.get("EGRI")),
      ndti_mean: safeNumber(feature.get("NDTI")),
      red_turbidity_proxy_mean: safeNumber(feature.get("RED_TURBIDITY_PROXY")),
      hossain_red_ntu_proxy_mean: safeNumber(feature.get("HOSSAIN_RED_NTU_PROXY"))
    }).select(exportFieldOrder);
  }).sort("distance_m");
}

function buildMonthConfigs(startYear, endYear, months) {
  var monthNames = {
    1: "January",
    2: "February",
    3: "March",
    4: "April",
    5: "May",
    6: "June",
    7: "July",
    8: "August",
    9: "September",
    10: "October",
    11: "November",
    12: "December"
  };

  var configs = [];
  for (var year = startYear; year <= endYear; year++) {
    months.forEach(function(month) {
      var monthLabel = (month < 10 ? "0" : "") + month;
      configs.push({
        label: year + "-" + monthLabel,
        title: monthNames[month] + " " + year,
        start: ee.Date.fromYMD(year, month, 1)
      });
    });
  }
  return configs;
}

// ================= GEOMETRY SETUP =================
var polygons = buildPolygons();
var aoi = polygons.geometry().bounds().buffer(aoiBufferMeters);
var aoiPixelCount = ee.Number(
  ee.Image.constant(1)
    .rename("aoi_px")
    .reduceRegion({
      reducer: ee.Reducer.count(),
      geometry: aoi,
      scale: analysisScaleSentinelMeters,
      bestEffort: true,
      maxPixels: reduceRegionMaxPixels,
      tileScale: summaryTileScale
    })
    .get("aoi_px")
);

print("AOI", aoi);
print("Polygon count", polygons.size());
print("Max downstream distance (m)", polygons.aggregate_max("distance_m"));
print("Minimum valid pixels", minValidPixels);
print("Max Sentinel scenes per month", maxSentinelScenesPerMonth);
print("AOI buffer (m)", aoiBufferMeters);
print("Summary tileScale", summaryTileScale);
print("reduceRegion maxPixels", reduceRegionMaxPixels);
print("NDSSI formula", "(Blue - NIR) / (Blue + NIR)");
print("EGRI formula", "Green / Red");
print("NDTI formula", "(Red - Green) / (Red + Green)");
print("Red turbidity proxy", "Red surface reflectance");
print("Hossain red-band proxy", "2677.2 * pow(red, 1.8562); exploratory, not locally calibrated");

// ================= SENTINEL-2 WORKFLOW =================
function scaleAndSelectS2(img) {
  var reflectance = img
    .select(["B2", "B3", "B4", "B8"], ["blue", "green", "red", "nir"])
    .multiply(0.0001);

  return reflectance
    .addBands(img.select(["QA60", "SCL"]))
    .copyProperties(img, ["system:time_start", "CLOUDY_PIXEL_PERCENTAGE"]);
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
      sensor_label: "Sentinel-2",
      source_label: "COPERNICUS/S2_SR_HARMONIZED",
      analysis_scale_m: analysisScaleSentinelMeters
    })
    .copyProperties(img, ["system:time_start", "CLOUDY_PIXEL_PERCENTAGE"]);
}

function addAoiClearFraction(img) {
  var validPx = ee.Number(
    ee.Image.constant(1)
      .updateMask(img.select("red").mask())
      .rename("valid_px")
      .reduceRegion({
        reducer: ee.Reducer.count(),
        geometry: aoi,
        scale: analysisScaleSentinelMeters,
        bestEffort: true,
        maxPixels: reduceRegionMaxPixels,
        tileScale: summaryTileScale
      })
      .get("valid_px")
  );

  var clearFraction = ee.Algorithms.If(aoiPixelCount.gt(0), validPx.divide(aoiPixelCount), 0);
  return img.set({
    aoi_valid_px: validPx,
    aoi_clear_fraction: clearFraction
  });
}

var sentinelBase = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  .filterBounds(aoi)
  .filterDate(
    ee.Date.fromYMD(sentinelStartYear, sentinelMonths[0], 1),
    ee.Date.fromYMD(sentinelEndYear, sentinelMonths[sentinelMonths.length - 1], 1).advance(1, "month")
  )
  .map(scaleAndSelectS2)
  .map(maskS2)
  .map(addCoreIndices)
  .map(addAoiClearFraction)
  .filter(ee.Filter.gte("aoi_clear_fraction", minSceneClearFractionOverAoi));

function buildSentinelMonthSummary(config) {
  var monthStart = config.start;
  var monthEnd = monthStart.advance(1, "month");
  var col = sentinelBase
    .filterDate(monthStart, monthEnd)
    .sort("aoi_clear_fraction", false)
    .limit(maxSentinelScenesPerMonth);
  var imageCount = col.size();
  var monthImage = ee.Image(ee.Algorithms.If(
    imageCount.gt(0),
    col.median(),
    emptyAnalysisImage()
  ));

  return summarizeImageToPolygons(monthImage, polygons, {
    export_group: "monthly_summary",
    dataset_label: "Sentinel-2 downstream monthly",
    sensor_label: "Sentinel-2",
    source_label: "COPERNICUS/S2_SR_HARMONIZED",
    scene_label: config.title,
    image_label: config.title,
    analysis_date: config.label,
    acquisition_date: null,
    analysis_year: ee.Number.parse(ee.String(config.label).slice(0, 4)),
    analysis_month: ee.Number.parse(ee.String(config.label).slice(5, 7)),
    month_key: config.label,
    month_label: config.title,
    window_start: monthStart.format("YYYY-MM-dd"),
    window_end_exclusive: monthEnd.format("YYYY-MM-dd"),
    image_count: imageCount,
    per_image_valid_count: imageCount,
    analysis_scale_m: analysisScaleSentinelMeters
  });
}

var sentinelStats = ee.FeatureCollection([]);
if (exportSentinelCsv || exportCombinedCsv) {
  var sentinelConfigs = buildMonthConfigs(sentinelStartYear, sentinelEndYear, sentinelMonths);
  var sentinelRaw = ee.FeatureCollection(sentinelConfigs.map(buildSentinelMonthSummary)).flatten();
  sentinelStats = selectExportFields(addRelativeFields(sentinelRaw, "month_key"));
  print("Sentinel export row count", sentinelStats.size());
}

// ================= PLANET / RAPIDEYE WORKFLOW =================
function makePlanetImage(scene) {
  var sceneDate = ee.Date(scene.acquisitionDate);
  var srImages = scene.srAssetIds.map(function(assetId) {
    return ee.Image(assetId)
      .select([0, 1, 2, 3, 4], ["blue", "green", "red", "red_edge", "nir"])
      .multiply(0.0001);
  });
  var udmMasks = scene.udmAssetIds.map(function(assetId) {
    return ee.Image(assetId).select([0], ["UDM"]).eq(0);
  });

  var srMosaic = ee.ImageCollection.fromImages(srImages).mosaic();
  var clearMask = ee.Image(ee.Algorithms.If(
    requirePlanetUdmClearMask,
    ee.ImageCollection.fromImages(udmMasks).mosaic(),
    ee.Image.constant(1)
  ));

  var validMask = clearMask
    .and(srMosaic.select("blue").gt(0))
    .and(srMosaic.select("green").gt(0))
    .and(srMosaic.select("red").gt(0))
    .and(srMosaic.select("nir").gt(0));

  return addCoreIndices(srMosaic.updateMask(validMask))
    .set({
      "system:time_start": sceneDate.millis(),
      acquisition_date: sceneDate.format("YYYY-MM-dd"),
      analysis_date: sceneDate.format("YYYY-MM-dd"),
      analysis_year: sceneDate.get("year"),
      analysis_month: sceneDate.get("month"),
      month_key: sceneDate.format("YYYY-MM"),
      month_label: sceneDate.format("YYYY-MM"),
      scene_label: scene.label,
      image_label: scene.label,
      sensor_label: scene.sensor,
      source_label: "Planet Labs uploaded asset",
      analysis_scale_m: analysisScalePlanetMeters,
      sr_asset_ids: scene.srAssetIds.join(","),
      udm_asset_ids: scene.udmAssetIds.join(",")
    });
}

function buildPlanetSceneSummary(scene) {
  var sceneImage = makePlanetImage(scene);
  var sceneDate = ee.Date(scene.acquisitionDate);
  return summarizeImageToPolygons(sceneImage, polygons, {
    export_group: "planet_scene",
    dataset_label: "Planet Labs downstream scene",
    sensor_label: scene.sensor,
    source_label: "Planet Labs uploaded asset",
    scene_label: scene.label,
    image_label: scene.label,
    analysis_date: sceneDate.format("YYYY-MM-dd"),
    acquisition_date: sceneDate.format("YYYY-MM-dd"),
    analysis_year: sceneDate.get("year"),
    analysis_month: sceneDate.get("month"),
    month_key: sceneDate.format("YYYY-MM"),
    month_label: sceneDate.format("YYYY-MM"),
    window_start: sceneDate.format("YYYY-MM-dd"),
    window_end_exclusive: sceneDate.advance(1, "day").format("YYYY-MM-dd"),
    image_count: 1,
    per_image_valid_count: 1,
    analysis_scale_m: analysisScalePlanetMeters
  });
}

var planetStats = ee.FeatureCollection([]);
if ((exportPlanetCsv || exportCombinedCsv) && includePlanetScenes) {
  var configuredPlanetScenes = filterConfiguredPlanetScenes(planetScenes);
  print("Planet scene config count", planetScenes.length);
  print("Configured Planet scene count", configuredPlanetScenes.length);
  print("Planet scenes", configuredPlanetScenes);
  if (configuredPlanetScenes.length > 0) {
    planetStats = selectExportFields(
      addRelativeFields(
        ee.FeatureCollection(configuredPlanetScenes.map(buildPlanetSceneSummary)).flatten(),
        "scene_label"
      )
    );
    print("Planet export row count", planetStats.size());
  } else {
    print("Planet export skipped: no fully configured scene assets found.");
  }
}

// ================= EXPORTS =================
if (exportSentinelCsv) {
  Export.table.toDrive({
    collection: sentinelStats,
    description: "downstream_recent_febmar_2017_2026_metrics_with_ndti",
    folder: exportFolder,
    fileNamePrefix: "downstream_recent_febmar_2017_2026_metrics_with_ndti",
    fileFormat: "CSV"
  });
}

if (exportPlanetCsv) {
  Export.table.toDrive({
    collection: planetStats,
    description: "downstream_planet_scene_metrics_with_ndti",
    folder: exportFolder,
    fileNamePrefix: "downstream_planet_scene_metrics_with_ndti",
    fileFormat: "CSV"
  });
}

if (exportCombinedCsv) {
  Export.table.toDrive({
    collection: sentinelStats.merge(planetStats),
    description: "downstream_planet_s2_combined_metrics_with_ndti",
    folder: exportFolder,
    fileNamePrefix: "downstream_planet_s2_combined_metrics_with_ndti",
    fileFormat: "CSV"
  });
}

print("Unified downstream export script ready.");
print("Sentinel export enabled", exportSentinelCsv);
print("Planet export enabled", exportPlanetCsv);
print("Combined export enabled", exportCombinedCsv);
