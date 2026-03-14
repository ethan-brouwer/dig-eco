/*
  FILE: gee/turbidity/san_sebastian_impact_pool_image_turbidity.js
  PURPOSE: Build lightweight monthly Sentinel-2 water-quality summaries for the impact_pool and upstream_control polygons.

  GEE IMPORTS REQUIRED
  - impact_pool
  - upstream_control
*/

// ================= USER SETTINGS =================
var siteName = "San Sebastian impact pool";
var analysisStartYear = 2017;
var cloudMax = 60;
var analysisScaleMeters = 10;
var useWetSeasonFilter = false;
var wetSeasonMonths = [5, 6, 7, 8, 9, 10];
var ndwiThreshold = -0.05;
var mndwiThreshold = -0.05;
var nirWaterMax = 0.15;
var swirWaterMax = 0.12;
var minWaterPixels = 3;
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
print("NDWI threshold", ndwiThreshold);
print("MNDWI threshold", mndwiThreshold);
print("NDMI formula", "(Green - NIR) / (Green + NIR)");
print("NDSSI formula", "(Blue - NIR) / (Blue + NIR)");
print("EGRI formula", "Green / Red");
print("NIR max", nirWaterMax);
print("SWIR1 max", swirWaterMax);
print("Min water pixels", minWaterPixels);
print("Comparison areas", comparisonAreas);

// ================= S2 HELPERS =================
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
  return img.select(["B2", "B3", "B4", "B8", "B11", "QA60", "SCL"], [
    "B2", "B3", "B4", "B8", "B11", "QA60", "SCL"
  ]).addBands(
    img.select(["B2", "B3", "B4", "B8", "B11"])
      .multiply(0.0001),
    null,
    true
  );
}

function addMonth(img) {
  return img.set("month", ee.Date(img.get("system:time_start")).get("month"));
}

function addIndices(img) {
  var ndwi = img.normalizedDifference(["B3", "B8"]).rename("NDWI");
  var mndwi = img.normalizedDifference(["B3", "B11"]).rename("MNDWI");
  var ndmi = img.normalizedDifference(["B3", "B8"]).rename("NDMI");
  var nsmi = img.expression(
    "(RED + GREEN - BLUE) / (RED + GREEN + BLUE)",
    {
      RED: img.select("B4"),
      GREEN: img.select("B3"),
      BLUE: img.select("B2")
    }
  ).rename("NSMI");
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
  return img.addBands([ndwi, mndwi, ndmi, nsmi, ndssi, egri]);
}

function waterMask(img) {
  return img.select("NDWI").gt(ndwiThreshold)
    .or(img.select("MNDWI").gt(mndwiThreshold))
    .or(
      img.select("B8").lt(nirWaterMax)
        .and(img.select("B11").lt(swirWaterMax))
    )
    .selfMask()
    .rename("WATER_MASK");
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

// ================= MONTHLY AGGREGATION =================
var monthlyStarts = monthStartList(
  ee.Date(startDate.format("YYYY-MM-01")),
  ee.Date(endDate.format("YYYY-MM-01")).advance(1, "month")
);

var monthlyStats = ee.FeatureCollection(monthlyStarts.map(function(mStart) {
  mStart = ee.Date(mStart);
  var mEnd = mStart.advance(1, "month");
  var col = s2.filterDate(mStart, mEnd);
  var imageCount = col.size();

  var monthImg = ee.Image(ee.Algorithms.If(
    imageCount.gt(0),
    col.median().clip(aoi),
    ee.Image.constant([0, 0, 0, 0, 0, 0]).rename(["NDWI", "MNDWI", "NDMI", "NSMI", "NDSSI", "EGRI"]).selfMask()
  ));

  var monthlyFeatures = comparisonAreas.map(function(area) {
    area = ee.Feature(area);
    var areaGeom = area.geometry();
    var reachId = area.get("reach_id");
    var reachType = area.get("reach_type");

    var emptyFeature = ee.Feature(null, {
      "system:time_start": mStart.millis(),
      date: mStart.format("YYYY-MM"),
      year: ee.Number.parse(mStart.format("YYYY")),
      year_label: mStart.format("YYYY"),
      month_num: ee.Number.parse(mStart.format("M")),
      month_label: mStart.format("MMM"),
      reach_id: reachId,
      reach_type: reachType,
      image_count: imageCount,
      water_px: 0,
      qa_flag: "no_images",
      ndwi_mean: noDataValue,
      mndwi_mean: noDataValue,
      ndmi_mean: noDataValue,
      nsmi_mean: noDataValue,
      ndssi_mean: noDataValue,
      egri_mean: noDataValue
    });

    return ee.Feature(ee.Algorithms.If(imageCount.gt(0), (function() {
      var water = waterMask(monthImg).clip(areaGeom);

      var waterPx = safeNumber(water.reduceRegion({
        reducer: ee.Reducer.count(),
        geometry: areaGeom,
        scale: analysisScaleMeters,
        bestEffort: true,
        maxPixels: 1e8
      }).get("WATER_MASK"));

      var stats = monthImg
        .select(["NDWI", "MNDWI", "NDMI", "NSMI", "NDSSI", "EGRI"])
        .updateMask(water)
        .reduceRegion({
          reducer: ee.Reducer.mean(),
          geometry: areaGeom,
          scale: analysisScaleMeters,
          bestEffort: true,
          maxPixels: 1e8
        });

      var isValid = waterPx.gte(minWaterPixels);

      return ee.Feature(null, {
        "system:time_start": mStart.millis(),
        date: mStart.format("YYYY-MM"),
        year: ee.Number.parse(mStart.format("YYYY")),
        year_label: mStart.format("YYYY"),
        month_num: ee.Number.parse(mStart.format("M")),
        month_label: mStart.format("MMM"),
        reach_id: reachId,
        reach_type: reachType,
        image_count: imageCount,
        water_px: waterPx,
        qa_flag: ee.Algorithms.If(isValid, "monthly_average", "low_water_px"),
        ndwi_mean: ee.Algorithms.If(isValid, safeNumber(stats.get("NDWI")), noDataValue),
        mndwi_mean: ee.Algorithms.If(isValid, safeNumber(stats.get("MNDWI")), noDataValue),
        ndmi_mean: ee.Algorithms.If(isValid, safeNumber(stats.get("NDMI")), noDataValue),
        nsmi_mean: ee.Algorithms.If(isValid, safeNumber(stats.get("NSMI")), noDataValue),
        ndssi_mean: ee.Algorithms.If(isValid, safeNumber(stats.get("NDSSI")), noDataValue),
        egri_mean: ee.Algorithms.If(isValid, safeNumber(stats.get("EGRI")), noDataValue)
      });
    })(), emptyFeature));
  });

  return monthlyFeatures;
})).flatten().sort("system:time_start");

print("Monthly turbidity table", monthlyStats);

var validMonthlyStats = monthlyStats.filter(ee.Filter.eq("qa_flag", "monthly_average"));
print("Valid month count", validMonthlyStats.size());

var highNsmiMonths = validMonthlyStats.sort("nsmi_mean", false).limit(20);
var lowNdwiMonths = validMonthlyStats.sort("ndwi_mean", true).limit(20);
var highNdssiMonths = validMonthlyStats.sort("ndssi_mean", false).limit(20);

print("Highest NSMI months", highNsmiMonths);
print("Lowest NDWI months", lowNdwiMonths);
print("Highest NDSSI months", highNdssiMonths);

// ================= CHARTS =================
var monthlyChartFc = monthlyStats.map(function(f) {
  return f.set({
    ndwi_mean_chart: cleanChartValue(f, "ndwi_mean"),
    mndwi_mean_chart: cleanChartValue(f, "mndwi_mean"),
    ndmi_mean_chart: cleanChartValue(f, "ndmi_mean"),
    nsmi_mean_chart: cleanChartValue(f, "nsmi_mean"),
    ndssi_mean_chart: cleanChartValue(f, "ndssi_mean"),
    egri_mean_chart: cleanChartValue(f, "egri_mean"),
    water_px_chart: ee.Algorithms.If(ee.Number(f.get("water_px")).gt(0), f.get("water_px"), null)
  });
});

var monthlyWaterChart = ui.Chart.feature.groups(
  monthlyChartFc.filter(ee.Filter.notNull(["ndwi_mean_chart"])),
  "date",
  "ndwi_mean_chart",
  "reach_id"
).setChartType("LineChart")
  .setOptions({
    title: "Monthly NDWI Comparison",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "NDWI"},
    lineWidth: 2,
    pointSize: 3,
    colors: ["#E53935", "#26A69A"]
  });
print(monthlyWaterChart);

var monthlyMndwiChart = ui.Chart.feature.groups(
  monthlyChartFc.filter(ee.Filter.notNull(["mndwi_mean_chart"])),
  "date",
  "mndwi_mean_chart",
  "reach_id"
).setChartType("LineChart")
  .setOptions({
    title: "Monthly MNDWI Comparison",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "MNDWI"},
    lineWidth: 2,
    pointSize: 3,
    colors: ["#E53935", "#26A69A"]
  });
print(monthlyMndwiChart);

var monthlyNdmiChart = ui.Chart.feature.groups(
  monthlyChartFc.filter(ee.Filter.notNull(["ndmi_mean_chart"])),
  "date",
  "ndmi_mean_chart",
  "reach_id"
).setChartType("LineChart")
  .setOptions({
    title: "Monthly NDMI Comparison",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "NDMI"},
    lineWidth: 2,
    pointSize: 3,
    colors: ["#E53935", "#26A69A"]
  });
print(monthlyNdmiChart);

var monthlyNsmiChart = ui.Chart.feature.groups(
  monthlyChartFc.filter(ee.Filter.notNull(["nsmi_mean_chart"])),
  "date",
  "nsmi_mean_chart",
  "reach_id"
).setChartType("LineChart")
  .setOptions({
    title: "Monthly NSMI Comparison",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "NSMI"},
    lineWidth: 2,
    pointSize: 3,
    colors: ["#E53935", "#26A69A"]
  });
print(monthlyNsmiChart);

var monthlyNdssiChart = ui.Chart.feature.groups(
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
print(monthlyNdssiChart);

var monthlyEgriChart = ui.Chart.feature.groups(
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
print(monthlyEgriChart);

function makeYearOverlayChart(fc, propertyName, title, axisTitle) {
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
      pointSize: 3
    });
}

var impactChartFc = monthlyChartFc.filter(ee.Filter.eq("reach_id", "impact_pool"));
var controlChartFc = monthlyChartFc.filter(ee.Filter.eq("reach_id", "upstream_control"));

print(makeYearOverlayChart(impactChartFc, "ndwi_mean_chart", "Impact Pool NDWI Monthly Seasonality by Year", "NDWI monthly mean"));
print(makeYearOverlayChart(controlChartFc, "ndwi_mean_chart", "Upstream Control NDWI Monthly Seasonality by Year", "NDWI monthly mean"));
print(makeYearOverlayChart(impactChartFc, "nsmi_mean_chart", "Impact Pool NSMI Monthly Seasonality by Year", "NSMI monthly mean"));
print(makeYearOverlayChart(controlChartFc, "nsmi_mean_chart", "Upstream Control NSMI Monthly Seasonality by Year", "NSMI monthly mean"));
print(makeYearOverlayChart(impactChartFc, "ndssi_mean_chart", "Impact Pool NDSSI Monthly Seasonality by Year", "NDSSI monthly mean"));
print(makeYearOverlayChart(controlChartFc, "ndssi_mean_chart", "Upstream Control NDSSI Monthly Seasonality by Year", "NDSSI monthly mean"));

var waterPxChart = ui.Chart.feature.groups(
  monthlyChartFc.filter(ee.Filter.notNull(["water_px_chart"])),
  "date",
  "water_px_chart",
  "reach_id"
).setChartType("ColumnChart")
  .setOptions({
    title: "Monthly Water Pixels Used",
    hAxis: {title: "Month", slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "Water pixels at 10 m"},
    colors: ["#E53935", "#26A69A"]
  });
print(waterPxChart);

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
  latestPreview.select("NSMI"),
  {min: -0.3, max: 0.3, palette: ["#2166AC", "#F7F7F7", "#B2182B"]},
  "Latest monthly NSMI",
  false
);
