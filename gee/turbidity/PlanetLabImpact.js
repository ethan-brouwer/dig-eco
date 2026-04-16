/*
  FILE: gee/turbidity/PlanetLabImpact.js
  PURPOSE: Add uploaded Planet Labs RapidEye tiles to the existing impact
  pool EGRI and NDSSI workflow while keeping the same imported GEE polygons.

  GEE IMPORTS REQUIRED
  - impact_pool
  - upstream_control

  ASSET SETUP
  1. Upload each RapidEye surface reflectance GeoTIFF as an Earth Engine image.

  2. Upload each matching RapidEye UDM GeoTIFF as an Earth Engine image.

  3. Replace the srAssetId and udmAssetId placeholders in planetScenes below.

  Current local scenes:
  - 2011-02-27 RapidEye-5
    SR:  turb/ImpactPool(extra)_reorthotile_analytic_sr/REOrthoTile/1645310_2011-02-27_RE5_3A_Analytic_SR_clip.tif
    UDM: turb/ImpactPool(extra)_reorthotile_analytic_sr/REOrthoTile/1645310_2011-02-27_RE5_3A_udm_clip.tif

  - 2012-02-14 RapidEye-4
    SR:  turb/ImpactPool_2/REOrthoTile/1645310_2012-02-14_RE4_3A_Analytic_SR_clip.tif
    UDM: turb/ImpactPool_2/REOrthoTile/1645310_2012-02-14_RE4_3A_udm_clip.tif

  RapidEye band order used here:
  B1 = blue, B2 = green, B3 = red, B4 = red edge, B5 = nir
*/

// ================= USER SETTINGS =================
var analysisStartYear = 2011;
var cloudMax = 60;
var sentinelScaleMeters = 10;
var planetScaleMeters = 5;
var comparisonMonth = 2;
var useWetSeasonFilter = false;
var wetSeasonMonths = [5, 6, 7, 8, 9, 10];
var minValidPixels = 3;
var noDataValue = -9999;

var includePlanetLabsTiles = true;
var planetScenes = [
  {
    label: "RapidEye-5 2011-02-27",
    sensor: "RapidEye-5",
    acquisitionDate: "2011-02-27",
    srAssetId: "projects/YOUR_PROJECT/assets/ImpactPool_extra_RE_SR_20110227",
    udmAssetId: "projects/YOUR_PROJECT/assets/ImpactPool_extra_RE_UDM_20110227"
  },
  {
    label: "RapidEye-4 2012-02-14",
    sensor: "RapidEye-4",
    acquisitionDate: "2012-02-14",
    srAssetId: "projects/YOUR_PROJECT/assets/ImpactPool_2_RE_SR_20120214",
    udmAssetId: "projects/YOUR_PROJECT/assets/ImpactPool_2_RE_UDM_20120214"
  }
];

// ================= GEOMETRY HELPERS =================
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

function cleanChartValue(f, propertyName) {
  var value = f.get(propertyName);
  return ee.Algorithms.If(ee.Algorithms.IsEqual(value, noDataValue), null, value);
}

function monthStartList(start, end) {
  var monthCount = ee.Number(end.difference(start, "month")).floor();
  return ee.List.sequence(0, monthCount.subtract(1)).map(function(m) {
    return start.advance(ee.Number(m), "month");
  });
}

function makeYearOverlayChart(fc, propertyName, title, axisTitle, colors) {
  return ui.Chart.feature.groups(
    fc.filter(ee.Filter.notNull([propertyName])),
    "month_num",
    propertyName,
    "year_label"
  ).setChartType("LineChart")
    .setOptions({
      title: title,
      hAxis: {
        title: "Month of year",
        ticks: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
      },
      vAxis: {title: axisTitle},
      lineWidth: 2,
      pointSize: 4,
      colors: colors
    });
}

var impactPoolGeom = toGeometry(impact_pool);
var upstreamControlGeom = toGeometry(upstream_control);
var aoi = impactPoolGeom.union(upstreamControlGeom, 1);

var comparisonAreas = ee.FeatureCollection([
  ee.Feature(impactPoolGeom, {
    reach_id: "impact_pool",
    reach_type: "impact"
  }),
  ee.Feature(upstreamControlGeom, {
    reach_id: "upstream_control",
    reach_type: "control"
  })
]);

// ================= MAP SETUP =================
Map.setOptions("SATELLITE");
Map.centerObject(impactPoolGeom, 16);
Map.addLayer(impactPoolGeom, {color: "FFFF00"}, "impact_pool", true);
Map.addLayer(upstreamControlGeom, {color: "26A69A"}, "upstream_control", true);

// ================= DATE WINDOW =================
var endDate = ee.Date(Date.now());
var startDate = ee.Date.fromYMD(analysisStartYear, 1, 1);

print("Analysis window", startDate.format("YYYY-MM-dd"), endDate.format("YYYY-MM-dd"));
print("Collection filter geometry", "impact_pool + upstream_control");
print("Sentinel-2 scale (m)", sentinelScaleMeters);
print("Planet Labs scale (m)", planetScaleMeters);
print("Comparison month", comparisonMonth);
print("Wet-season month filter", useWetSeasonFilter ? wetSeasonMonths : "OFF");
print("Min valid pixels", minValidPixels);
print("NDSSI formula", "(Blue - NIR) / (Blue + NIR)");
print("EGRI formula", "Green / Red");
print("Comparison areas", comparisonAreas);

// ================= INDEX HELPERS =================
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

  return img.addBands([ndssi, egri]);
}

function addMonth(img) {
  return img.set({
    month: ee.Date(img.get("system:time_start")).get("month"),
    year: ee.Date(img.get("system:time_start")).get("year")
  });
}

// ================= SENTINEL-2 HELPERS =================
function prepS2(img) {
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

  return img
    .updateMask(bad.not())
    .select(["blue", "green", "red", "nir"])
    .set({
      sensor: "Sentinel-2",
      source: "COPERNICUS/S2_SR_HARMONIZED",
      analysis_scale_m: sentinelScaleMeters
    })
    .copyProperties(img, ["system:time_start", "CLOUDY_PIXEL_PERCENTAGE"]);
}

var s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  .filterBounds(aoi)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", cloudMax))
  .map(prepS2)
  .map(maskS2)
  .map(addIndices)
  .map(addMonth);

if (useWetSeasonFilter) {
  s2 = s2.filter(ee.Filter.inList("month", wetSeasonMonths));
}

// ================= PLANET LABS RAPIDEYE HELPERS =================
function makePlanetLabsImage(scene) {
  var sceneDate = ee.Date(scene.acquisitionDate);
  var sr = ee.Image(scene.srAssetId)
    .select([0, 1, 2, 3, 4], ["blue", "green", "red", "red_edge", "nir"])
    .multiply(0.0001);

  var udm = ee.Image(scene.udmAssetId).select([0], ["UDM"]);

  var validMask = udm.eq(0)
    .and(sr.select("blue").gt(0))
    .and(sr.select("green").gt(0))
    .and(sr.select("red").gt(0))
    .and(sr.select("nir").gt(0));

  return addIndices(sr.updateMask(validMask))
    .set({
      "system:time_start": sceneDate.millis(),
      date: sceneDate.format("YYYY-MM-dd"),
      month: sceneDate.get("month"),
      year: sceneDate.get("year"),
      scene_label: scene.label,
      sensor: scene.sensor,
      source: "Planet Labs uploaded asset",
      sr_asset_id: scene.srAssetId,
      udm_asset_id: scene.udmAssetId,
      analysis_scale_m: planetScaleMeters
    });
}

var planetCollection = ee.ImageCollection([]);

if (includePlanetLabsTiles) {
  planetCollection = ee.ImageCollection(planetScenes.map(makePlanetLabsImage));

  print("Planet Labs scene config", planetScenes);
  print("Planet Labs image count", planetCollection.size());
  print("Planet Labs observations", planetCollection);

  planetScenes.forEach(function(scene, index) {
    var sceneImage = makePlanetLabsImage(scene).clip(aoi);
    var shown = index === planetScenes.length - 1;

    Map.addLayer(
      sceneImage.select(["red", "green", "blue"]),
      {min: 0.02, max: 0.35, gamma: 1.2},
      "Planet Labs RGB " + scene.label,
      shown
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
}

var analysisImages = s2.merge(planetCollection).sort("system:time_start");

print("Sentinel-2 image count", s2.size());
print("Combined image count", analysisImages.size());

// ================= MONTHLY AREA SUMMARIES =================
var monthlyStarts = monthStartList(
  ee.Date(startDate.format("YYYY-MM-01")),
  ee.Date(endDate.format("YYYY-MM-01")).advance(1, "month")
);

var monthlyStats = ee.FeatureCollection(monthlyStarts.map(function(mStart) {
  mStart = ee.Date(mStart);
  var mEnd = mStart.advance(1, "month");
  var col = analysisImages.filterDate(mStart, mEnd);
  var imageCount = col.size();
  var dateLabel = mStart.format("YYYY-MM");

  var monthlyFeatures = comparisonAreas.map(function(area) {
    area = ee.Feature(area);
    var areaGeom = area.geometry();

    var emptyFeature = ee.Feature(null, {
      "system:time_start": mStart.millis(),
      date: dateLabel,
      year: ee.Number.parse(mStart.format("YYYY")),
      year_label: mStart.format("YYYY"),
      month_num: ee.Number.parse(mStart.format("M")),
      reach_id: area.get("reach_id"),
      reach_type: area.get("reach_type"),
      image_count: imageCount,
      sensors: "none",
      analysis_scale_m: sentinelScaleMeters,
      valid_px: 0,
      qa_flag: "no_images",
      ndssi_mean: noDataValue,
      egri_mean: noDataValue
    });

    return ee.Feature(ee.Algorithms.If(imageCount.gt(0), (function() {
      var monthImg = ee.Image(col.median()).clip(areaGeom);
      var validMask = monthImg.select("red").mask().clip(areaGeom);
      var sensorNames = col.aggregate_array("sensor").distinct().sort().join(", ");
      var scale = ee.Number(col.aggregate_min("analysis_scale_m"));

      var validPx = safeNumber(validMask.reduceRegion({
        reducer: ee.Reducer.count(),
        geometry: areaGeom,
        scale: scale,
        bestEffort: true,
        maxPixels: 1e8
      }).get("red"));

      var stats = monthImg
        .select(["NDSSI", "EGRI"])
        .updateMask(validMask)
        .reduceRegion({
          reducer: ee.Reducer.mean(),
          geometry: areaGeom,
          scale: scale,
          bestEffort: true,
          maxPixels: 1e8
        });

      var hasValidData = validPx.gte(minValidPixels);

      return ee.Feature(null, {
        "system:time_start": mStart.millis(),
        date: dateLabel,
        year: ee.Number.parse(mStart.format("YYYY")),
        year_label: mStart.format("YYYY"),
        month_num: ee.Number.parse(mStart.format("M")),
        reach_id: area.get("reach_id"),
        reach_type: area.get("reach_type"),
        image_count: imageCount,
        sensors: sensorNames,
        analysis_scale_m: scale,
        valid_px: validPx,
        qa_flag: ee.Algorithms.If(hasValidData, "valid", "low_valid_px"),
        ndssi_mean: ee.Algorithms.If(hasValidData, safeNumber(stats.get("NDSSI")), noDataValue),
        egri_mean: ee.Algorithms.If(hasValidData, safeNumber(stats.get("EGRI")), noDataValue)
      });
    })(), emptyFeature));
  });

  return monthlyFeatures;
})).flatten().sort("system:time_start");

var validMonthlyStats = monthlyStats.filter(ee.Filter.eq("qa_flag", "valid"));

print("Monthly comparison table", monthlyStats);
print("Valid month count", validMonthlyStats.size());

// ================= IMPACT-CONTROL DELTAS =================
var deltaStats = ee.FeatureCollection(monthlyStarts.map(function(mStart) {
  mStart = ee.Date(mStart);
  var dateLabel = mStart.format("YYYY-MM");
  var monthNum = ee.Number.parse(mStart.format("M"));
  var yearLabel = mStart.format("YYYY");

  var impactFc = validMonthlyStats
    .filter(ee.Filter.eq("date", dateLabel))
    .filter(ee.Filter.eq("reach_id", "impact_pool"));
  var controlFc = validMonthlyStats
    .filter(ee.Filter.eq("date", dateLabel))
    .filter(ee.Filter.eq("reach_id", "upstream_control"));

  var hasBoth = impactFc.size().gt(0).and(controlFc.size().gt(0));
  var impact = ee.Feature(ee.Algorithms.If(
    hasBoth,
    impactFc.first(),
    ee.Feature(null, {egri_mean: noDataValue, ndssi_mean: noDataValue})
  ));
  var control = ee.Feature(ee.Algorithms.If(
    hasBoth,
    controlFc.first(),
    ee.Feature(null, {egri_mean: noDataValue, ndssi_mean: noDataValue})
  ));

  return ee.Feature(null, {
    "system:time_start": mStart.millis(),
    date: dateLabel,
    year: ee.Number.parse(mStart.format("YYYY")),
    year_label: yearLabel,
    month_num: monthNum,
    qa_flag: ee.Algorithms.If(hasBoth, "paired_valid", "missing_pair"),
    sensors: ee.Algorithms.If(hasBoth, impact.get("sensors"), "none"),
    egri_impact: ee.Algorithms.If(hasBoth, impact.get("egri_mean"), noDataValue),
    egri_control: ee.Algorithms.If(hasBoth, control.get("egri_mean"), noDataValue),
    egri_delta: ee.Algorithms.If(
      hasBoth,
      ee.Number(impact.get("egri_mean")).subtract(ee.Number(control.get("egri_mean"))),
      noDataValue
    ),
    ndssi_impact: ee.Algorithms.If(hasBoth, impact.get("ndssi_mean"), noDataValue),
    ndssi_control: ee.Algorithms.If(hasBoth, control.get("ndssi_mean"), noDataValue),
    ndssi_delta: ee.Algorithms.If(
      hasBoth,
      ee.Number(impact.get("ndssi_mean")).subtract(ee.Number(control.get("ndssi_mean"))),
      noDataValue
    )
  });
})).sort("system:time_start");

var pairedDeltaStats = deltaStats.filter(ee.Filter.eq("qa_flag", "paired_valid"));

print("Monthly impact-minus-control deltas", deltaStats);
print("Largest positive EGRI delta months", pairedDeltaStats.sort("egri_delta", false).limit(15));
print("Largest positive NDSSI delta months", pairedDeltaStats.sort("ndssi_delta", false).limit(15));

// ================= CHARTS =================
var monthlyChartFc = monthlyStats.map(function(f) {
  return f.set({
    egri_mean_chart: cleanChartValue(f, "egri_mean"),
    ndssi_mean_chart: cleanChartValue(f, "ndssi_mean"),
    valid_px_chart: ee.Algorithms.If(ee.Number(f.get("valid_px")).gt(0), f.get("valid_px"), null)
  });
});

var deltaChartFc = deltaStats.map(function(f) {
  return f.set({
    egri_delta_chart: cleanChartValue(f, "egri_delta"),
    ndssi_delta_chart: cleanChartValue(f, "ndssi_delta")
  });
});

var februaryStats = monthlyChartFc
  .filter(ee.Filter.eq("month_num", comparisonMonth))
  .filter(ee.Filter.eq("qa_flag", "valid"));

var februaryDeltas = deltaChartFc
  .filter(ee.Filter.eq("month_num", comparisonMonth))
  .filter(ee.Filter.eq("qa_flag", "paired_valid"));

print("February comparison table", februaryStats);
print("February impact-minus-control deltas", februaryDeltas);

var februaryEgriChart = ui.Chart.feature.groups(
  februaryStats.filter(ee.Filter.notNull(["egri_mean_chart"])),
  "year_label",
  "egri_mean_chart",
  "reach_id"
).setChartType("ColumnChart")
  .setOptions({
    title: "February EGRI by Year",
    hAxis: {title: "Year"},
    vAxis: {title: "Mean EGRI"},
    bar: {groupWidth: "75%"},
    colors: ["#E53935", "#26A69A"],
    legend: {position: "top"}
  });
print(februaryEgriChart);

var februaryNdssiChart = ui.Chart.feature.groups(
  februaryStats.filter(ee.Filter.notNull(["ndssi_mean_chart"])),
  "year_label",
  "ndssi_mean_chart",
  "reach_id"
).setChartType("ColumnChart")
  .setOptions({
    title: "February NDSSI by Year",
    hAxis: {title: "Year"},
    vAxis: {title: "Mean NDSSI"},
    bar: {groupWidth: "75%"},
    colors: ["#E53935", "#26A69A"],
    legend: {position: "top"}
  });
print(februaryNdssiChart);

var februaryEgriDeltaChart = ui.Chart.feature.byFeature(
  februaryDeltas.filter(ee.Filter.notNull(["egri_delta_chart"])),
  "year_label",
  ["egri_delta_chart"]
).setChartType("ColumnChart")
  .setOptions({
    title: "February EGRI Delta: impact_pool - upstream_control",
    hAxis: {title: "Year"},
    vAxis: {title: "EGRI delta"},
    bar: {groupWidth: "75%"},
    colors: ["#8E24AA"]
  });
print(februaryEgriDeltaChart);

var februaryNdssiDeltaChart = ui.Chart.feature.byFeature(
  februaryDeltas.filter(ee.Filter.notNull(["ndssi_delta_chart"])),
  "year_label",
  ["ndssi_delta_chart"]
).setChartType("ColumnChart")
  .setOptions({
    title: "February NDSSI Delta: impact_pool - upstream_control",
    hAxis: {title: "Year"},
    vAxis: {title: "NDSSI delta"},
    bar: {groupWidth: "75%"},
    colors: ["#3949AB"]
  });
print(februaryNdssiDeltaChart);

var februaryValidPxChart = ui.Chart.feature.groups(
  februaryStats.filter(ee.Filter.notNull(["valid_px_chart"])),
  "year_label",
  "valid_px_chart",
  "reach_id"
).setChartType("ColumnChart")
  .setOptions({
    title: "February Valid Pixels Used",
    hAxis: {title: "Year"},
    vAxis: {title: "Valid pixels"},
    bar: {groupWidth: "75%"},
    colors: ["#E53935", "#26A69A"],
    legend: {position: "top"}
  });
print(februaryValidPxChart);

// ================= QUICKLOOKS =================
var latestMonth = ee.Date(endDate.format("YYYY-MM-01"));
var latestMonthCol = s2.filterDate(latestMonth.advance(-1, "month"), latestMonth);
var latestPreview = ee.Image(ee.Algorithms.If(
  latestMonthCol.size().gt(0),
  latestMonthCol.median(),
  s2.sort("system:time_start", false).first()
)).clip(aoi);

Map.addLayer(
  latestPreview.select(["red", "green", "blue"]),
  {min: 0.02, max: 0.3},
  "Latest Sentinel-2 monthly RGB",
  false
);
Map.addLayer(
  latestPreview.select("EGRI"),
  {min: 0.5, max: 2.0, palette: ["#7F0000", "#FDD49E", "#238B45"]},
  "Latest Sentinel-2 monthly EGRI",
  false
);
Map.addLayer(
  latestPreview.select("NDSSI"),
  {min: -0.8, max: 0.4, palette: ["#54278F", "#2B8CBE", "#F7FCF0"]},
  "Latest Sentinel-2 monthly NDSSI",
  false
);

// ================= EXPORTS =================
Export.table.toDrive({
  collection: monthlyStats,
  description: "planet_lab_impact_monthly_egri_ndssi",
  fileFormat: "CSV"
});

Export.table.toDrive({
  collection: deltaStats,
  description: "planet_lab_impact_monthly_egri_ndssi_deltas",
  fileFormat: "CSV"
});

Export.table.toDrive({
  collection: februaryStats,
  description: "planet_lab_impact_february_egri_ndssi",
  fileFormat: "CSV"
});

Export.table.toDrive({
  collection: februaryDeltas,
  description: "planet_lab_impact_february_egri_ndssi_deltas",
  fileFormat: "CSV"
});
