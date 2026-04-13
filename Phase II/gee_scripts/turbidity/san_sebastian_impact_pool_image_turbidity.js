/*
  FILE: gee/turbidity/san_sebastian_impact_pool_image_turbidity.js
  PURPOSE: Compare lightweight monthly EGRI and NDSSI summaries for impact_pool and upstream_control.

  GEE IMPORTS REQUIRED
  - impact_pool
  - upstream_control
*/

// ================= USER SETTINGS =================
var analysisStartYear = 2017;
var cloudMax = 60;
var analysisScaleMeters = 10;
var useWetSeasonFilter = false;
var wetSeasonMonths = [5, 6, 7, 8, 9, 10];
var minValidPixels = 3;
var noDataValue = -9999;

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
      pointSize: 3,
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
print("Scale (m)", analysisScaleMeters);
print("Wet-season month filter", useWetSeasonFilter ? wetSeasonMonths : "OFF");
print("Min valid pixels", minValidPixels);
print("NDSSI formula", "(Blue - NIR) / (Blue + NIR)");
print("EGRI formula", "Green / Red");
print("Comparison areas", comparisonAreas);

// ================= S2 HELPERS =================
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

function addMonth(img) {
  return img.set("month", ee.Date(img.get("system:time_start")).get("month"));
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

var s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  .filterBounds(aoi)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", cloudMax))
  .map(scaleAndSelectS2)
  .map(maskS2)
  .map(addMonth)
  .map(addIndices);

if (useWetSeasonFilter) {
  s2 = s2.filter(ee.Filter.inList("month", wetSeasonMonths));
}

print("Raw image count", s2.size());

// ================= MONTHLY AREA SUMMARIES =================
var monthlyStarts = monthStartList(
  ee.Date(startDate.format("YYYY-MM-01")),
  ee.Date(endDate.format("YYYY-MM-01")).advance(1, "month")
);

var monthlyStats = ee.FeatureCollection(monthlyStarts.map(function(mStart) {
  mStart = ee.Date(mStart);
  var mEnd = mStart.advance(1, "month");
  var col = s2.filterDate(mStart, mEnd);
  var imageCount = col.size();

  var monthlyFeatures = comparisonAreas.map(function(area) {
    area = ee.Feature(area);
    var areaGeom = area.geometry();

    var emptyFeature = ee.Feature(null, {
      "system:time_start": mStart.millis(),
      date: mStart.format("YYYY-MM"),
      year: ee.Number.parse(mStart.format("YYYY")),
      year_label: mStart.format("YYYY"),
      month_num: ee.Number.parse(mStart.format("M")),
      reach_id: area.get("reach_id"),
      reach_type: area.get("reach_type"),
      image_count: imageCount,
      valid_px: 0,
      qa_flag: "no_images",
      ndssi_mean: noDataValue,
      egri_mean: noDataValue
    });

    return ee.Feature(ee.Algorithms.If(imageCount.gt(0), (function() {
      var monthImg = ee.Image(col.median()).clip(areaGeom);
      var validMask = monthImg.select("B4").mask().clip(areaGeom);

      var validPx = safeNumber(validMask.reduceRegion({
        reducer: ee.Reducer.count(),
        geometry: areaGeom,
        scale: analysisScaleMeters,
        bestEffort: true,
        maxPixels: 1e8
      }).get("B4"));

      var stats = monthImg
        .select(["NDSSI", "EGRI"])
        .updateMask(validMask)
        .reduceRegion({
          reducer: ee.Reducer.mean(),
          geometry: areaGeom,
          scale: analysisScaleMeters,
          bestEffort: true,
          maxPixels: 1e8
        });

      var hasValidData = validPx.gte(minValidPixels);

      return ee.Feature(null, {
        "system:time_start": mStart.millis(),
        date: mStart.format("YYYY-MM"),
        year: ee.Number.parse(mStart.format("YYYY")),
        year_label: mStart.format("YYYY"),
        month_num: ee.Number.parse(mStart.format("M")),
        reach_id: area.get("reach_id"),
        reach_type: area.get("reach_type"),
        image_count: imageCount,
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

var egriComparisonChart = ui.Chart.feature.groups(
  monthlyChartFc.filter(ee.Filter.notNull(["egri_mean_chart"])),
  "date",
  "egri_mean_chart",
  "reach_id"
).setChartType("LineChart")
  .setOptions({
    title: "Monthly EGRI Comparison",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "EGRI"},
    lineWidth: 2,
    pointSize: 3,
    colors: ["#E53935", "#26A69A"]
  });
print(egriComparisonChart);

var ndssiComparisonChart = ui.Chart.feature.groups(
  monthlyChartFc.filter(ee.Filter.notNull(["ndssi_mean_chart"])),
  "date",
  "ndssi_mean_chart",
  "reach_id"
).setChartType("LineChart")
  .setOptions({
    title: "Monthly NDSSI Comparison",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "NDSSI"},
    lineWidth: 2,
    pointSize: 3,
    colors: ["#E53935", "#26A69A"]
  });
print(ndssiComparisonChart);

var egriDeltaChart = ui.Chart.feature.byFeature(
  deltaChartFc.filter(ee.Filter.notNull(["egri_delta_chart"])),
  "date",
  ["egri_delta_chart"]
).setChartType("LineChart")
  .setOptions({
    title: "Monthly EGRI Delta: impact_pool - upstream_control",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "EGRI delta"},
    lineWidth: 2,
    pointSize: 3,
    colors: ["#8E24AA"]
  });
print(egriDeltaChart);

var ndssiDeltaChart = ui.Chart.feature.byFeature(
  deltaChartFc.filter(ee.Filter.notNull(["ndssi_delta_chart"])),
  "date",
  ["ndssi_delta_chart"]
).setChartType("LineChart")
  .setOptions({
    title: "Monthly NDSSI Delta: impact_pool - upstream_control",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "NDSSI delta"},
    lineWidth: 2,
    pointSize: 3,
    colors: ["#3949AB"]
  });
print(ndssiDeltaChart);

var impactChartFc = monthlyChartFc.filter(ee.Filter.eq("reach_id", "impact_pool"));
var controlChartFc = monthlyChartFc.filter(ee.Filter.eq("reach_id", "upstream_control"));

print(makeYearOverlayChart(
  impactChartFc,
  "egri_mean_chart",
  "Impact Pool EGRI Monthly Seasonality by Year",
  "EGRI monthly mean"
));
print(makeYearOverlayChart(
  controlChartFc,
  "egri_mean_chart",
  "Upstream Control EGRI Monthly Seasonality by Year",
  "EGRI monthly mean"
));
print(makeYearOverlayChart(
  impactChartFc,
  "ndssi_mean_chart",
  "Impact Pool NDSSI Monthly Seasonality by Year",
  "NDSSI monthly mean"
));
print(makeYearOverlayChart(
  controlChartFc,
  "ndssi_mean_chart",
  "Upstream Control NDSSI Monthly Seasonality by Year",
  "NDSSI monthly mean"
));
print(makeYearOverlayChart(
  deltaChartFc,
  "egri_delta_chart",
  "EGRI Delta Monthly Seasonality by Year",
  "EGRI delta"
));
print(makeYearOverlayChart(
  deltaChartFc,
  "ndssi_delta_chart",
  "NDSSI Delta Monthly Seasonality by Year",
  "NDSSI delta"
));

var validPxChart = ui.Chart.feature.groups(
  monthlyChartFc.filter(ee.Filter.notNull(["valid_px_chart"])),
  "date",
  "valid_px_chart",
  "reach_id"
).setChartType("ColumnChart")
  .setOptions({
    title: "Monthly Valid Pixels Used",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "Valid pixels at 10 m"},
    colors: ["#E53935", "#26A69A"]
  });
print(validPxChart);

// ================= QUICKLOOKS =================
var latestMonth = ee.Date(endDate.format("YYYY-MM-01"));
var latestMonthCol = s2.filterDate(latestMonth.advance(-1, "month"), latestMonth);
var latestPreview = ee.Image(ee.Algorithms.If(
  latestMonthCol.size().gt(0),
  latestMonthCol.median(),
  s2.sort("system:time_start", false).first()
)).clip(aoi);

Map.addLayer(
  latestPreview.select(["B4", "B3", "B2"]),
  {min: 0.02, max: 0.3},
  "Latest monthly RGB",
  false
);
Map.addLayer(
  latestPreview.select("EGRI"),
  {min: 0.5, max: 2.0, palette: ["#7F0000", "#FDD49E", "#238B45"]},
  "Latest monthly EGRI",
  false
);
