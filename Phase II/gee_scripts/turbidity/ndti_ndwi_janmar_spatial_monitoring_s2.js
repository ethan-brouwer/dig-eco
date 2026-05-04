/*
  FILE: Phase II/gee_scripts/turbidity/ndti_ndwi_janmar_spatial_monitoring_s2.js
  PURPOSE: Build January-March annual Sentinel-2 composites for spatial
           monitoring of water-only NDTI hotspots. The workflow keeps scenes
           with low tile-level cloudiness, masks cloudy/shadowed pixels,
           computes a median Jan-Mar composite for each year, masks to water
           using NDWI > 0.1, and exports pixel-level NDTI/NDWI values for
           graphing and hotspot analysis in Jupyter.

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
var analysisEndYear = new Date().getFullYear();
var compositeStartMonth = 1; // January
var compositeDurationMonths = 3; // January through March
var tileCloudMax = 30; // Relaxed tile-level prefilter; pixel QA mask still removes clouds
var compositeScaleMeters = 10;
var ndwiWaterThreshold = 0.0; // More permissive for narrow or mixed river pixels
var hotspotPercentile = 90;
var exportFolder = "EarthEngine";
var exportPixelTable = false;
var exportAnnualSummaryTable = false;
var exportHotspotRasters = false;
var showPreviewCharts = true;
var chartPreviewYear = 2024;

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
    .copyProperties(img, ["system:time_start", "CLOUDY_PIXEL_PERCENTAGE"]);
}

function addWaterAndSedimentIndices(img) {
  var ndwi = img.normalizedDifference(["green", "nir"]).rename("NDWI");
  var ndti = img.normalizedDifference(["red", "green"]).rename("NDTI");

  return img.addBands([ndwi, ndti]);
}

function makeEmptyComposite() {
  return ee.Image.constant([0, 0, 0, 0, 0, 0])
    .rename(["blue", "green", "red", "nir", "NDWI", "NDTI"])
    .selfMask();
}

function buildAnnualConfig(year) {
  var start = ee.Date.fromYMD(year, compositeStartMonth, 1);
  return {
    year: year,
    season_label: year + " Jan-Mar",
    start: start,
    end: start.advance(compositeDurationMonths, "month")
  };
}

function buildAnnualComposite(config) {
  var s2 = s2Base
    .filterDate(config.start, config.end)
    .sort("system:time_start");

  var imageCount = s2.size();
  var composite = ee.Image(ee.Algorithms.If(
    imageCount.gt(0),
    s2.median(),
    makeEmptyComposite()
  ));

  var waterMask = composite.select("NDWI").gt(ndwiWaterThreshold);
  var waterComposite = composite.updateMask(waterMask);

  var waterPixelCount = ee.Number(ee.Algorithms.If(
    imageCount.gt(0),
    ee.Image.constant(1)
      .updateMask(waterComposite.select("NDTI").mask())
      .rename("water_px")
      .reduceRegion({
        reducer: ee.Reducer.count(),
        geometry: aoi,
        scale: compositeScaleMeters,
        bestEffort: true,
        maxPixels: 1e8,
        tileScale: 2
      })
      .get("water_px"),
    0
  ));

  var percentileStats = ee.Dictionary(ee.Algorithms.If(
    imageCount.gt(0).and(waterPixelCount.gt(0)),
    waterComposite.select("NDTI").reduceRegion({
      reducer: ee.Reducer.percentile([50, 75, hotspotPercentile, 95]),
      geometry: aoi,
      scale: compositeScaleMeters,
      bestEffort: true,
      maxPixels: 1e8,
      tileScale: 2
    }),
    ee.Dictionary({})
  ));

  var hotspotThreshold = ee.Algorithms.If(
    waterPixelCount.gt(0),
    percentileStats.get("NDTI_p" + hotspotPercentile),
    null
  );
  var hotspotMask = ee.Image(ee.Algorithms.If(
    ee.Algorithms.IsEqual(hotspotThreshold, null),
    ee.Image(0).selfMask(),
    waterComposite.select("NDTI").gte(ee.Number(hotspotThreshold)).selfMask()
  )).rename("NDTI_hotspot");

  var sampleImage = waterComposite.select(["NDWI", "NDTI", "red", "green", "blue", "nir"])
    .addBands(ee.Image.constant(config.year).rename("analysis_year"))
    .addBands(ee.Image.pixelLonLat());

  var pixelSamples = ee.FeatureCollection(ee.Algorithms.If(
    imageCount.gt(0).and(waterPixelCount.gt(0)),
    sampleImage.sample({
      region: aoi,
      scale: compositeScaleMeters,
      geometries: true,
      tileScale: 2
    }).map(function(feature) {
      feature = ee.Feature(feature);
      return feature.set({
        composite_year: config.year,
        season_label: config.season_label,
        composite_start: config.start.format("YYYY-MM-dd"),
        composite_end_exclusive: config.end.format("YYYY-MM-dd"),
        image_count: imageCount,
        ndwi_water_threshold: ndwiWaterThreshold,
        hotspot_percentile: hotspotPercentile,
        hotspot_ndti_threshold: hotspotThreshold
      });
    }),
    ee.FeatureCollection([])
  ));

  var annualSummary = ee.Feature(null, {
    composite_year: config.year,
    composite_year_num: ee.Number(config.year),
    season_label: config.season_label,
    composite_start: config.start.format("YYYY-MM-dd"),
    composite_end_exclusive: config.end.format("YYYY-MM-dd"),
    image_count: imageCount,
    ndwi_water_threshold: ndwiWaterThreshold,
    hotspot_percentile: hotspotPercentile,
    hotspot_ndti_threshold: hotspotThreshold,
    water_pixel_count: waterPixelCount,
    ndti_p50: ee.Algorithms.If(waterPixelCount.gt(0), percentileStats.get("NDTI_p50"), null),
    ndti_p75: ee.Algorithms.If(waterPixelCount.gt(0), percentileStats.get("NDTI_p75"), null),
    ndti_p90: ee.Algorithms.If(waterPixelCount.gt(0), percentileStats.get("NDTI_p" + hotspotPercentile), null),
    ndti_p95: ee.Algorithms.If(waterPixelCount.gt(0), percentileStats.get("NDTI_p95"), null)
  });

  var polygonWaterDiagnostics = ee.FeatureCollection(ee.Algorithms.If(
    imageCount.gt(0),
    composite.select("NDWI").reduceRegions({
      collection: polygons,
      reducer: ee.Reducer.mean().combine({
        reducer2: ee.Reducer.max(),
        sharedInputs: true
      }),
      scale: compositeScaleMeters,
      tileScale: 2
    }),
    polygons.map(function(feature) {
      return ee.Feature(feature).set({
        mean: null,
        max: null
      });
    })
  )).map(function(feature) {
    feature = ee.Feature(feature);
    var geom = feature.geometry();
    var totalPolygonPx = ee.Number(
      ee.Image.constant(1)
        .rename("polygon_px")
        .reduceRegion({
          reducer: ee.Reducer.count(),
          geometry: geom,
          scale: compositeScaleMeters,
          bestEffort: true,
          maxPixels: 1e8,
          tileScale: 2
        })
        .get("polygon_px")
    );

    var qualifyingWaterPx = ee.Number(ee.Algorithms.If(
      imageCount.gt(0),
      ee.Image.constant(1)
        .updateMask(waterMask)
        .rename("water_px")
        .reduceRegion({
          reducer: ee.Reducer.count(),
          geometry: geom,
          scale: compositeScaleMeters,
          bestEffort: true,
          maxPixels: 1e8,
          tileScale: 2
        })
        .get("water_px"),
      0
    ));

    return feature.set({
      composite_year: config.year,
      composite_year_num: ee.Number(config.year),
      season_label: config.season_label,
      image_count: imageCount,
      ndwi_water_threshold: ndwiWaterThreshold,
      polygon_px_total: totalPolygonPx,
      polygon_water_px: qualifyingWaterPx,
      polygon_water_fraction: ee.Algorithms.If(
        totalPolygonPx.gt(0),
        qualifyingWaterPx.divide(totalPolygonPx),
        null
      ),
      polygon_qualifies_as_water: qualifyingWaterPx.gt(0),
      polygon_ndwi_mean: feature.get("mean"),
      polygon_ndwi_max: feature.get("max")
    });
  });

  return {
    year: config.year,
    imageCount: imageCount,
    composite: composite,
    waterComposite: waterComposite,
    hotspotMask: hotspotMask,
    pixelSamples: pixelSamples,
    annualSummary: annualSummary,
    polygonWaterDiagnostics: polygonWaterDiagnostics
  };
}

// ================= INPUTS =================
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
  ee.Feature(toGeometry(Poly5000m), {polygon_id: "Poly5000m", distance_m: 5000})
]);

var aoi = polygons.geometry().bounds().buffer(250);
var polygonUnion = polygons.geometry();
var outsidePolygons = aoi.difference(polygonUnion, 1);

print("Analysis years", analysisStartYear + " to " + analysisEndYear);
print("Composite start month", compositeStartMonth);
print("Composite duration (months)", compositeDurationMonths);
print("Tile cloud max (%)", tileCloudMax);
print("Composite scale (m)", compositeScaleMeters);
print("NDWI water threshold", ndwiWaterThreshold);
print("NDTI hotspot percentile", hotspotPercentile);
print("Show preview charts", showPreviewCharts);
print("Chart preview year", chartPreviewYear);
print("NDWI formula", "(GREEN - NIR) / (GREEN + NIR)");
print("NDTI formula", "(RED - GREEN) / (RED + GREEN)");

var s2Base = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  .filterBounds(aoi)
  .filterDate(
    ee.Date.fromYMD(analysisStartYear, compositeStartMonth, 1),
    ee.Date.fromYMD(analysisEndYear + 1, compositeStartMonth, 1).advance(compositeDurationMonths, "month")
  )
  .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", tileCloudMax))
  .map(scaleAndSelectS2)
  .map(maskS2)
  .map(addWaterAndSedimentIndices);

// ================= OUTPUTS =================
Map.setOptions("SATELLITE");
Map.centerObject(aoi, 14);
Map.addLayer(polygons, {color: "FFB300"}, "Monitoring polygons", true);

var annualOutputs = [];
var sampleCollections = [];
var annualSummaryFeatures = [];
var polygonDiagnosticCollections = [];

for (var year = analysisStartYear; year <= analysisEndYear; year++) {
  var annualOutput = buildAnnualComposite(buildAnnualConfig(year));
  annualOutputs.push(annualOutput);
  sampleCollections.push(annualOutput.pixelSamples);
  annualSummaryFeatures.push(annualOutput.annualSummary);
  polygonDiagnosticCollections.push(annualOutput.polygonWaterDiagnostics);

  Map.addLayer(
    annualOutput.composite.select(["red", "green", "blue"]),
    {min: 0.02, max: 0.3},
    year + " Jan-Mar RGB",
    year === analysisStartYear
  );
  Map.addLayer(
    annualOutput.waterComposite.select("NDWI"),
    {min: -0.2, max: 0.6, palette: ["#8c510a", "#f6e8c3", "#01665e"]},
    year + " Jan-Mar NDWI water mask",
    false
  );
  Map.addLayer(
    annualOutput.waterComposite.select("NDTI"),
    {min: -0.2, max: 0.4, palette: ["#313695", "#ffffbf", "#a50026"]},
    year + " Jan-Mar NDTI",
    false
  );
  Map.addLayer(
    annualOutput.hotspotMask,
    {palette: ["#ff00ff"]},
    year + " Jan-Mar NDTI hotspots",
    false
  );

  print(year + " image count", annualOutput.imageCount);
  print(year + " polygon water diagnostics", annualOutput.polygonWaterDiagnostics);
}

var allPixelSamples = ee.FeatureCollection([]);
sampleCollections.forEach(function(fc) {
  allPixelSamples = allPixelSamples.merge(fc);
});

var annualSummaryTable = ee.FeatureCollection(annualSummaryFeatures);
var polygonWaterDiagnosticsTable = ee.FeatureCollection([]);
polygonDiagnosticCollections.forEach(function(fc) {
  polygonWaterDiagnosticsTable = polygonWaterDiagnosticsTable.merge(fc);
});

print("Annual summary table", annualSummaryTable);
print("Polygon water diagnostics table", polygonWaterDiagnosticsTable);
print("Pixel sample row count", allPixelSamples.size());
print("Pixel sample preview", allPixelSamples.limit(10));

if (showPreviewCharts) {
  var annualSummaryChart = ui.Chart.feature.byFeature({
    features: annualSummaryTable.sort("composite_year"),
    xProperty: "composite_year_num",
    yProperties: ["ndti_p50", "ndti_p75", "ndti_p90", "ndti_p95"]
  }).setChartType("LineChart").setOptions({
    title: "Annual Jan-Mar Water-Only NDTI Percentiles",
    hAxis: {title: "Year"},
    vAxis: {title: "NDTI"},
    lineWidth: 2,
    pointSize: 5,
    series: {
      0: {color: "#2c7fb8"},
      1: {color: "#41b6c4"},
      2: {color: "#f03b20"},
      3: {color: "#bd0026"}
    }
  });
  print(annualSummaryChart);

  var waterPixelChart = ui.Chart.feature.byFeature({
    features: annualSummaryTable.sort("composite_year"),
    xProperty: "composite_year_num",
    yProperties: ["water_pixel_count", "image_count"]
  }).setChartType("ColumnChart").setOptions({
    title: "Annual Jan-Mar Water Pixel Count and Scene Count",
    hAxis: {title: "Year"},
    vAxis: {title: "Count"},
    series: {
      0: {targetAxisIndex: 0, color: "#1b9e77"},
      1: {targetAxisIndex: 1, color: "#7570b3"}
    },
    vAxes: {
      0: {title: "Water Pixels"},
      1: {title: "Scenes"}
    }
  });
  print(waterPixelChart);

  var polygonWaterChart = ui.Chart.feature.groups({
    features: polygonWaterDiagnosticsTable,
    xProperty: "composite_year_num",
    yProperty: "polygon_water_fraction",
    seriesProperty: "polygon_id"
  }).setChartType("LineChart").setOptions({
    title: "Jan-Mar Water Fraction by Polygon",
    hAxis: {title: "Year"},
    vAxis: {title: "Water-qualified fraction of polygon"},
    lineWidth: 2,
    pointSize: 4
  });
  print(polygonWaterChart);

  var polygonNdwiChart = ui.Chart.feature.groups({
    features: polygonWaterDiagnosticsTable,
    xProperty: "composite_year_num",
    yProperty: "polygon_ndwi_mean",
    seriesProperty: "polygon_id"
  }).setChartType("LineChart").setOptions({
    title: "Jan-Mar Mean NDWI by Polygon",
    hAxis: {title: "Year"},
    vAxis: {title: "Mean NDWI"},
    lineWidth: 2,
    pointSize: 4
  });
  print(polygonNdwiChart);

  var previewOutput = null;
  annualOutputs.forEach(function(output) {
    if (output.year === chartPreviewYear) {
      previewOutput = output;
    }
  });

  if (previewOutput !== null) {
    var ndtiHistogram = ui.Chart.image.histogram({
      image: previewOutput.waterComposite.select("NDTI"),
      region: aoi,
      scale: compositeScaleMeters,
      maxPixels: 1e8
    }).setOptions({
      title: chartPreviewYear + " Jan-Mar Water-Only NDTI Histogram",
      hAxis: {title: "NDTI"},
      vAxis: {title: "Pixel count"},
      colors: ["#d95f02"]
    });
    print(ndtiHistogram);

    var ndwiPreMaskHistogram = ui.Chart.image.histogram({
      image: previewOutput.composite.select("NDWI"),
      region: aoi,
      scale: compositeScaleMeters,
      maxPixels: 1e8
    }).setOptions({
      title: chartPreviewYear + " Jan-Mar NDWI Histogram Before Water Mask",
      hAxis: {title: "NDWI"},
      vAxis: {title: "Pixel count"},
      colors: ["#6a3d9a"]
    });
    print(ndwiPreMaskHistogram);

    var ndwiInsidePolygonsHistogram = ui.Chart.image.histogram({
      image: previewOutput.composite.select("NDWI"),
      region: polygonUnion,
      scale: compositeScaleMeters,
      maxPixels: 1e8
    }).setOptions({
      title: chartPreviewYear + " Jan-Mar NDWI Histogram Inside Monitoring Polygons",
      hAxis: {title: "NDWI"},
      vAxis: {title: "Pixel count"},
      colors: ["#33a02c"]
    });
    print(ndwiInsidePolygonsHistogram);

    var ndwiOutsidePolygonsHistogram = ui.Chart.image.histogram({
      image: previewOutput.composite.select("NDWI"),
      region: outsidePolygons,
      scale: compositeScaleMeters,
      maxPixels: 1e8
    }).setOptions({
      title: chartPreviewYear + " Jan-Mar NDWI Histogram Outside Monitoring Polygons",
      hAxis: {title: "NDWI"},
      vAxis: {title: "Pixel count"},
      colors: ["#e31a1c"]
    });
    print(ndwiOutsidePolygonsHistogram);

    var ndwiHistogram = ui.Chart.image.histogram({
      image: previewOutput.waterComposite.select("NDWI"),
      region: aoi,
      scale: compositeScaleMeters,
      maxPixels: 1e8
    }).setOptions({
      title: chartPreviewYear + " Jan-Mar Water-Only NDWI Histogram",
      hAxis: {title: "NDWI"},
      vAxis: {title: "Pixel count"},
      colors: ["#1f78b4"]
    });
    print(ndwiHistogram);
  } else {
    print("Chart preview year not found in annual outputs:", chartPreviewYear);
  }
}

if (exportPixelTable) {
  Export.table.toDrive({
    collection: allPixelSamples,
    description: "ndti_ndwi_janmar_water_pixels_s2",
    folder: exportFolder,
    fileNamePrefix: "ndti_ndwi_janmar_water_pixels_s2",
    fileFormat: "CSV"
  });
}

if (exportAnnualSummaryTable) {
  Export.table.toDrive({
    collection: annualSummaryTable,
    description: "ndti_ndwi_janmar_annual_summary_s2",
    folder: exportFolder,
    fileNamePrefix: "ndti_ndwi_janmar_annual_summary_s2",
    fileFormat: "CSV"
  });
}

if (exportHotspotRasters) {
  annualOutputs.forEach(function(output) {
    Export.image.toDrive({
      image: output.hotspotMask.clip(aoi),
      description: "ndti_hotspots_janmar_" + output.year,
      folder: exportFolder,
      fileNamePrefix: "ndti_hotspots_janmar_" + output.year,
      region: aoi,
      scale: compositeScaleMeters,
      maxPixels: 1e9
    });
  });
}

print("Script ready. It builds annual Jan-Mar water-only composites, computes NDWI and NDTI, and previews hotspot behavior in GEE.");
