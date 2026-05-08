/*
  FILE: Phase I/gee_scripts/turbidity/san_sebastian_downstream_red_blue_monthly_fast.js
  PURPOSE: Fast downstream 0-10 km export for monthly 3-image composites.

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
  - Poly6000m
  - Poly7000m
  - Poly8000m
  - Poly9000m
  - Poly10000m

  OUTPUT
  - One CSV export with monthly median blue/red values per polygon.
  - One in-GEE annual NDTI preview chart with one line per year.

  NOTES
  - This script keeps the load light by:
    1. using Sentinel-2 only,
    2. limiting each month to the 3 least-cloudy scenes,
    3. exporting only red/blue medians plus small metadata.
  - The preview NDTI here is based on RED and BLUE:
      (red_median - blue_median) / (red_median + blue_median)
*/

// ================= USER SETTINGS =================
var analysisStartYear = 2017;
var analysisEndYear = 2025; // Keep complete years by default.
var cloudMax = 60;
var imagesPerMonth = 3;
var analysisScaleMeters = 10;
var exportFolder = "GEE_exports";
var exportPrefix = "san_sebastian_downstream_red_blue_monthly_3img_2017_2025";

// Set to a subset like [2, 3] if you only want Feb/Mar.
var monthsToUse = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];

var noDataValue = -9999;

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
  var cloud = qa.bitwiseAnd(1 << 10).neq(0).or(qa.bitwiseAnd(1 << 11).neq(0));
  var shadow = scl.eq(3);
  var cirrus = scl.eq(10);
  var snow = scl.eq(11);
  var saturated = scl.eq(1);
  var bad = cloud.or(shadow).or(cirrus).or(snow).or(saturated);
  return img.updateMask(bad.not());
}

function scaleAndSelectS2(img) {
  return img.select(["B2", "B4"], ["blue", "red"])
    .multiply(0.0001)
    .copyProperties(img, ["system:time_start", "CLOUDY_PIXEL_PERCENTAGE"]);
}

function addDateParts(img) {
  var d = ee.Date(img.get("system:time_start"));
  return img.set({
    year: d.get("year"),
    month: d.get("month"),
    month_key: d.format("YYYY-MM")
  });
}

function safeNumber(x) {
  return ee.Number(ee.Algorithms.If(ee.Algorithms.IsEqual(x, null), noDataValue, x));
}

function monthStartList(startYear, endYear, months) {
  var starts = [];
  for (var y = startYear; y <= endYear; y++) {
    months.forEach(function(m) {
      starts.push(ee.Date.fromYMD(y, m, 1));
    });
  }
  return ee.List(starts);
}

function makeEmptyFeature(monthStart, feature, qaFlag) {
  feature = ee.Feature(feature);
  return ee.Feature(null, {
    date: monthStart.format("YYYY-MM"),
    year: ee.Number.parse(monthStart.format("YYYY")),
    month_num: ee.Number.parse(monthStart.format("M")),
    month_label: monthStart.format("MMM"),
    polygon_id: feature.get("polygon_id"),
    distance_m: feature.get("distance_m"),
    image_count_available: 0,
    image_count_used: 0,
    qa_flag: qaFlag,
    blue_median: noDataValue,
    red_median: noDataValue
  });
}

// ================= INPUT POLYGONS =================
var polygons = ee.FeatureCollection([
  ee.Feature(toGeometry(upstream_control), {polygon_id: "upstream_control", distance_m: -500}),
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
  ee.Feature(toGeometry(Poly5000m), {polygon_id: "Poly5000m", distance_m: 5000}),
  ee.Feature(toGeometry(Poly6000m), {polygon_id: "Poly6000m", distance_m: 6000}),
  ee.Feature(toGeometry(Poly7000m), {polygon_id: "Poly7000m", distance_m: 7000}),
  ee.Feature(toGeometry(Poly8000m), {polygon_id: "Poly8000m", distance_m: 8000}),
  ee.Feature(toGeometry(Poly9000m), {polygon_id: "Poly9000m", distance_m: 9000}),
  ee.Feature(toGeometry(Poly10000m), {polygon_id: "Poly10000m", distance_m: 10000})
]).sort("distance_m");

var aoi = polygons.geometry().bounds().buffer(250);

// ================= COLLECTION =================
var startDate = ee.Date.fromYMD(analysisStartYear, 1, 1);
var endDate = ee.Date.fromYMD(analysisEndYear + 1, 1, 1);

var s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  .filterBounds(aoi)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", cloudMax))
  .map(maskS2)
  .map(scaleAndSelectS2)
  .map(addDateParts);

print("Analysis window", startDate.format("YYYY-MM-dd"), endDate.advance(-1, "day").format("YYYY-MM-dd"));
print("Months included", monthsToUse);
print("Image count after broad filters", s2.size());

Map.setOptions("SATELLITE");
Map.centerObject(aoi, 14);
Map.addLayer(polygons, {color: "FFB300"}, "Downstream 0-10 km polygons", true);

// ================= MONTHLY EXPORT TABLE =================
var monthlyStarts = monthStartList(analysisStartYear, analysisEndYear, monthsToUse);

var monthlyLong = ee.FeatureCollection(monthlyStarts.map(function(monthStart) {
  monthStart = ee.Date(monthStart);
  var monthEnd = monthStart.advance(1, "month");

  var monthAll = s2.filterDate(monthStart, monthEnd);
  var monthTop = monthAll
    .sort("CLOUDY_PIXEL_PERCENTAGE")
    .limit(imagesPerMonth);

  var availableCount = monthAll.size();
  var usedCount = monthTop.size();

  var monthComposite = ee.Image(ee.Algorithms.If(
    usedCount.gt(0),
    monthTop.median().clip(aoi),
    ee.Image.constant([0, 0]).rename(["blue", "red"]).selfMask()
  ));

  return polygons.map(function(feature) {
    feature = ee.Feature(feature);
    var geom = feature.geometry();

    return ee.Feature(ee.Algorithms.If(usedCount.gt(0), (function() {
      var stats = monthComposite.reduceRegion({
        reducer: ee.Reducer.median(),
        geometry: geom,
        scale: analysisScaleMeters,
        maxPixels: 1e8,
        tileScale: 2
      });

      var blueMedian = safeNumber(stats.get("blue"));
      var redMedian = safeNumber(stats.get("red"));
      var hasData = ee.Algorithms.If(
        blueMedian.eq(noDataValue),
        false,
        ee.Algorithms.If(redMedian.eq(noDataValue), false, true)
      );

      return ee.Feature(null, {
        date: monthStart.format("YYYY-MM"),
        year: ee.Number.parse(monthStart.format("YYYY")),
        month_num: ee.Number.parse(monthStart.format("M")),
        month_label: monthStart.format("MMM"),
        polygon_id: feature.get("polygon_id"),
        distance_m: feature.get("distance_m"),
        image_count_available: availableCount,
        image_count_used: usedCount,
        qa_flag: ee.Algorithms.If(hasData, "ok", "no_valid_pixels"),
        blue_median: blueMedian,
        red_median: redMedian
      });
    })(), makeEmptyFeature(monthStart, feature, "no_images")));
  });
})).flatten();

print("Monthly export table", monthlyLong.limit(25));

// ================= PREVIEW CHART: ONE LINE PER YEAR =================
var monthlyForPreview = monthlyLong.map(function(f) {
  var blue = ee.Number(f.get("blue_median"));
  var red = ee.Number(f.get("red_median"));
  var valid = ee.Algorithms.If(
    blue.eq(noDataValue),
    false,
    ee.Algorithms.If(
      red.eq(noDataValue),
      false,
      ee.Algorithms.If(red.add(blue).eq(0), false, true)
    )
  );
  var ndtiPreview = ee.Algorithms.If(
    valid,
    red.subtract(blue).divide(red.add(blue)),
    null
  );
  return ee.Feature(f).set("ndti_rb_preview", ndtiPreview);
});

var yearList = ee.List.sequence(analysisStartYear, analysisEndYear);
var annualProfiles = ee.FeatureCollection(yearList.map(function(y) {
  y = ee.Number(y);
  return polygons.map(function(poly) {
    poly = ee.Feature(poly);
    var subset = monthlyForPreview
      .filter(ee.Filter.eq("year", y))
      .filter(ee.Filter.eq("polygon_id", poly.get("polygon_id")))
      .filter(ee.Filter.notNull(["ndti_rb_preview"]));

    var annualNdti = ee.Algorithms.If(
      subset.size().gt(0),
      ee.Array(subset.aggregate_array("ndti_rb_preview"))
        .reduce(ee.Reducer.median(), [0])
        .get([0]),
      null
    );

    return ee.Feature(null, {
      year: y,
      polygon_id: poly.get("polygon_id"),
      distance_m: poly.get("distance_m"),
      ndti_rb_year: annualNdti
    });
  });
})).flatten().filter(ee.Filter.notNull(["ndti_rb_year"]));

print(ui.Chart.feature.groups(annualProfiles, "distance_m", "ndti_rb_year", "year")
  .setChartType("LineChart")
  .setOptions({
    title: "Annual RED/BLUE NDTI Preview by Distance",
    hAxis: {title: "Distance downstream (m)"},
    vAxis: {title: "(Red - Blue) / (Red + Blue)"},
    lineWidth: 2,
    pointSize: 3
  }));

// ================= LIGHT MAP PREVIEW =================
var previewMonth = s2
  .filter(ee.Filter.calendarRange(analysisEndYear, analysisEndYear, "year"))
  .filter(ee.Filter.calendarRange(monthsToUse[0], monthsToUse[0], "month"))
  .sort("CLOUDY_PIXEL_PERCENTAGE")
  .limit(imagesPerMonth)
  .median()
  .clip(aoi);

Map.addLayer(previewMonth, {bands: ["red", "red", "blue"], min: 0.02, max: 0.30}, "Preview monthly composite", false);

// ================= EXPORT =================
Export.table.toDrive({
  collection: monthlyLong,
  description: exportPrefix,
  folder: exportFolder,
  fileNamePrefix: exportPrefix,
  fileFormat: "CSV",
  selectors: [
    "date",
    "year",
    "month_num",
    "month_label",
    "polygon_id",
    "distance_m",
    "image_count_available",
    "image_count_used",
    "qa_flag",
    "blue_median",
    "red_median"
  ]
});
