/*
  FILE: Phase II/gee_scripts/turbidity/downstream_recent_janfeb_2017_2026.js
  PURPOSE: Analyze January-February downstream optical water-quality proxies
           for every year from 2017 through 2026 using recent Sentinel-2
           imagery and the existing upstream/downstream polygons, with CSV
           export as the primary output for downstream analysis in Jupyter.

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
*/

// ================= USER SETTINGS =================
var analysisStartYear = 2017;
var analysisEndYear = 2026;
var analysisMonths = [1, 2]; // January and February
var cloudMax = 40;
var sentinelScaleMeters = 10;
var analysisScaleMeters = sentinelScaleMeters;
var upstreamReferenceDistanceMeters = -500;
var minValidPixels = 3;
var exportCsvTables = true;
var exportFolder = "EarthEngine";
var showPreviewCharts = false;

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
  var redMedian = ee.Algorithms.If(
    validPerImageStats.size().gt(0),
    ee.Array(validPerImageStats.aggregate_array("red_turbidity_proxy_mean"))
      .reduce(ee.Reducer.median(), [0])
      .get([0]),
    null
  );
  var hossainMedian = ee.Algorithms.If(
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
    has_valid_pixels: validPx.gte(minValidPixels),
    qa_flag: ee.Algorithms.If(imageCount.eq(0), "no_images",
      ee.Algorithms.If(validPx.gte(minValidPixels), "ok", "low_valid_px")),
    ndssi_mean: ndssiMedian,
    egri_mean: egriMedian,
    red_turbidity_proxy_mean: redMedian,
    hossain_red_ntu_proxy_mean: hossainMedian,
    per_image_valid_count: validPerImageStats.size()
  });
}

function normalizeByImpactPool(fc, groupProperty, valueProperty, outputProperty) {
  var groupValues = ee.List(fc.aggregate_array(groupProperty)).distinct();
  var normalized = ee.FeatureCollection(groupValues.map(function(groupValue) {
    var groupFc = fc.filter(ee.Filter.eq(groupProperty, groupValue));
    var baselineFeature = ee.Feature(groupFc.filter(ee.Filter.eq("polygon_id", "impact_pool")).first());
    var baselineRaw = baselineFeature.get(valueProperty);
    var baselineValue = ee.Number(baselineRaw);

    return groupFc.map(function(feature) {
      feature = ee.Feature(feature);
      var valueRaw = feature.get(valueProperty);
      var value = ee.Number(valueRaw);
      var normalizedValue = ee.Algorithms.If(
        ee.Algorithms.IsEqual(baselineRaw, null),
        null,
        ee.Algorithms.If(
          ee.Algorithms.IsEqual(valueRaw, null),
          null,
          ee.Algorithms.If(baselineValue.neq(0), value.divide(baselineValue), null)
        )
      );
      return feature.set(outputProperty, normalizedValue);
    });
  })).flatten();
  return normalized.filter(ee.Filter.notNull([outputProperty]));
}

function attachChartSeries(fc, groupProperty) {
  return fc.map(function(feature) {
    feature = ee.Feature(feature);
    var groupLabel = ee.String(feature.get(groupProperty));
    var chartSeries = ee.String(ee.Algorithms.If(
      ee.String(feature.get("polygon_id")).compareTo("upstream_control").eq(0),
      groupLabel.cat(" upstream"),
      groupLabel
    ));
    return feature.set("chart_series", chartSeries);
  });
}

function attachMatchedProperty(baseFc, lookupFc, groupProperty, lookupProperty, outputProperty) {
  return baseFc.map(function(feature) {
    feature = ee.Feature(feature);
    var groupValue = feature.get(groupProperty);
    var polygonId = feature.get("polygon_id");
    var match = ee.Feature(lookupFc
      .filter(ee.Filter.eq(groupProperty, groupValue))
      .filter(ee.Filter.eq("polygon_id", polygonId))
      .first());
    var matchedValue = ee.Algorithms.If(
      ee.Algorithms.IsEqual(match, null),
      null,
      match.get(lookupProperty)
    );
    return feature.set(outputProperty, matchedValue);
  });
}

function buildAnalysisConfigs(startYear, endYear, months) {
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

print("Sentinel-2 scale (m)", sentinelScaleMeters);
print("Analysis summary scale (m)", analysisScaleMeters);
print("Cloud max (%)", cloudMax);
print("Analysis years", analysisStartYear + " to " + analysisEndYear);
print("Analysis months", analysisMonths);
print("Upstream reference distance (m)", upstreamReferenceDistanceMeters);
print("Minimum valid pixels", minValidPixels);
print("Show preview charts", showPreviewCharts);
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
    .map(addIndices)
    .sort("system:time_start");

  var imageCount = s2.size();
  var monthImage = ee.Image(ee.Algorithms.If(
    imageCount.gt(0),
    s2.median().clip(aoi),
    makeEmptyAnalysisImage()
  ));

  var polygonStats = polygons.map(function(feature) {
    return summarizePolygon(feature, s2, monthImage, imageCount)
      .set({
        month_label: config.title,
        month_key: config.label,
        analysis_date: config.label
      });
  }).sort("distance_m");

  var validSignalStats = polygonStats.filter(ee.Filter.gte("valid_px", minValidPixels));

  return {
    title: config.title,
    label: config.label,
    monthStart: monthStart,
    monthEnd: monthEnd,
    imageCount: imageCount,
    monthImage: monthImage,
    polygonStats: polygonStats,
    validSignalStats: validSignalStats
  };
}

// ================= OUTPUTS =================
Map.setOptions("SATELLITE");
Map.centerObject(aoi, 15);
Map.addLayer(polygons, {color: "FFB300"}, "Downstream polygons", true);
Map.addLayer(toGeometry(upstream_control), {color: "26A69A"}, "upstream_control", true);
Map.addLayer(toGeometry(impact_pool), {color: "FFFF00"}, "impact_pool", true);

var analysisConfigs = buildAnalysisConfigs(analysisStartYear, analysisEndYear, analysisMonths);
print("Analysis windows", analysisConfigs);

var monthSummaries = analysisConfigs.map(buildMonthSummary);
var combinedMonthlyStats = ee.FeatureCollection(monthSummaries.map(function(summary) {
  return summary.polygonStats;
})).flatten();

var normalizedMonthlyEgriStats = normalizeByImpactPool(
  combinedMonthlyStats, "month_label", "egri_mean", "egri_rel_impact_pool"
);
var normalizedMonthlyNdssiStats = normalizeByImpactPool(
  combinedMonthlyStats, "month_label", "ndssi_mean", "ndssi_rel_impact_pool"
);
var normalizedMonthlyHossainStats = normalizeByImpactPool(
  combinedMonthlyStats, "month_label", "hossain_red_ntu_proxy_mean", "hossain_rel_impact_pool"
);
var normalizedMonthlyRedStats = normalizeByImpactPool(
  combinedMonthlyStats, "month_label", "red_turbidity_proxy_mean", "red_rel_impact_pool"
);

var monthlyExportStats = attachMatchedProperty(
  attachMatchedProperty(
    attachMatchedProperty(
      attachMatchedProperty(
        combinedMonthlyStats,
        normalizedMonthlyEgriStats,
        "month_label",
        "egri_rel_impact_pool",
        "egri_rel_impact_pool"
      ),
      normalizedMonthlyNdssiStats,
      "month_label",
      "ndssi_rel_impact_pool",
      "ndssi_rel_impact_pool"
    ),
    normalizedMonthlyHossainStats,
    "month_label",
    "hossain_rel_impact_pool",
    "hossain_rel_impact_pool"
  ),
  normalizedMonthlyRedStats,
  "month_label",
  "red_rel_impact_pool",
  "red_rel_impact_pool"
).map(function(feature) {
  feature = ee.Feature(feature);
  return feature.set({
    export_group: "monthly_janfeb_recent",
    analysis_year: ee.Number.parse(ee.String(feature.get("month_key")).slice(0, 4)),
    analysis_month: ee.Number.parse(ee.String(feature.get("month_key")).slice(5, 7)),
    polygon_role: ee.Algorithms.If(
      ee.String(feature.get("polygon_id")).compareTo("upstream_control").eq(0),
      "upstream_control",
      ee.Algorithms.If(
        ee.String(feature.get("polygon_id")).compareTo("impact_pool").eq(0),
        "impact_pool",
        "downstream"
      )
    ),
    sensor_label: "Sentinel-2",
    source_label: "COPERNICUS/S2_SR_HARMONIZED",
    analysis_scale_m: analysisScaleMeters
  });
});

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
  print(summary.title + " window", summary.monthStart.format("YYYY-MM-dd"), summary.monthEnd.format("YYYY-MM-dd"));
  print(summary.title + " image count", summary.imageCount);
  print(summary.title + " exportable polygon rows", summary.polygonStats.size());
});

print("Monthly export row count", monthlyExportStats.size());
print("Upstream rows in export", monthlyExportStats.filter(
  ee.Filter.eq("polygon_id", "upstream_control")
).size());

if (showPreviewCharts) {
  var chartMonthlyStats = attachChartSeries(
    combinedMonthlyStats.filter(ee.Filter.gte("valid_px", minValidPixels)),
    "month_label"
  );
  var chartNormalizedMonthlyEgriStats = attachChartSeries(normalizedMonthlyEgriStats, "month_label");
  var chartNormalizedMonthlyNdssiStats = attachChartSeries(normalizedMonthlyNdssiStats, "month_label");
  var chartNormalizedMonthlyHossainStats = attachChartSeries(normalizedMonthlyHossainStats, "month_label");

  print(ui.Chart.feature.groups(
    chartMonthlyStats, "distance_m", "egri_mean", "chart_series"
  ).setChartType("LineChart").setOptions({
    title: "January-February EGRI by Distance Downstream (2017-2026)",
    hAxis: {title: "Distance from impact pool (m)"},
    vAxis: {title: "Mean EGRI"},
    lineWidth: 2,
    pointSize: 4
  }));

  print(ui.Chart.feature.groups(
    chartMonthlyStats, "distance_m", "ndssi_mean", "chart_series"
  ).setChartType("LineChart").setOptions({
    title: "January-February NDSSI by Distance Downstream (2017-2026)",
    hAxis: {title: "Distance from impact pool (m)"},
    vAxis: {title: "Mean NDSSI"},
    lineWidth: 2,
    pointSize: 4
  }));

  print(ui.Chart.feature.groups(
    chartMonthlyStats, "distance_m", "hossain_red_ntu_proxy_mean", "chart_series"
  ).setChartType("LineChart").setOptions({
    title: "January-February Hossain Proxy by Distance Downstream (2017-2026)",
    hAxis: {title: "Distance from impact pool (m)"},
    vAxis: {title: "Exploratory NTU proxy"},
    lineWidth: 2,
    pointSize: 4
  }));

  print(ui.Chart.feature.groups(
    chartMonthlyStats, "distance_m", "red_turbidity_proxy_mean", "chart_series"
  ).setChartType("LineChart").setOptions({
    title: "January-February Red Reflectance Proxy by Distance Downstream (2017-2026)",
    hAxis: {title: "Distance from impact pool (m)"},
    vAxis: {title: "Mean red reflectance"},
    lineWidth: 2,
    pointSize: 4
  }));

  print(ui.Chart.feature.groups(
    chartNormalizedMonthlyEgriStats, "distance_m", "egri_rel_impact_pool", "chart_series"
  ).setChartType("LineChart").setOptions({
    title: "Relative EGRI by Distance (impact_pool = 1; 2017-2026 Jan-Feb)",
    hAxis: {title: "Distance from impact pool (m)"},
    vAxis: {title: "EGRI / impact_pool EGRI"},
    lineWidth: 2,
    pointSize: 4
  }));

  print(ui.Chart.feature.groups(
    chartNormalizedMonthlyNdssiStats, "distance_m", "ndssi_rel_impact_pool", "chart_series"
  ).setChartType("LineChart").setOptions({
    title: "Relative NDSSI by Distance (impact_pool = 1; 2017-2026 Jan-Feb)",
    hAxis: {title: "Distance from impact pool (m)"},
    vAxis: {title: "NDSSI / impact_pool NDSSI"},
    lineWidth: 2,
    pointSize: 4
  }));

  print(ui.Chart.feature.groups(
    chartNormalizedMonthlyHossainStats, "distance_m", "hossain_rel_impact_pool", "chart_series"
  ).setChartType("LineChart").setOptions({
    title: "Relative Hossain Proxy by Distance (impact_pool = 1; 2017-2026 Jan-Feb)",
    hAxis: {title: "Distance from impact pool (m)"},
    vAxis: {title: "Hossain proxy / impact_pool proxy"},
    lineWidth: 2,
    pointSize: 4
  }));
}

if (exportCsvTables) {
  Export.table.toDrive({
    collection: monthlyExportStats,
    description: "downstream_recent_janfeb_2017_2026_metrics",
    folder: exportFolder,
    fileNamePrefix: "downstream_recent_janfeb_2017_2026_metrics",
    fileFormat: "CSV"
  });
}

print("Script ready. It analyzes every January and February from 2017 through 2026 using recent Sentinel-2 imagery that passes the cloud and QA filters.");
