/*
  FILE: gee/groundwork/mrds_mine_disturbance_timeseries.js
  PURPOSE: Draft workflow to monitor disturbance around MRDS sites in El Salvador
           using Landsat annual composites and NDVI class thresholds.
  INPUTS: MRDS point asset, Landsat C2 L2 SR, GAUL boundary, DEM.
  OUTPUTS: Long and wide CSV tables for 1 km and 2 km buffers.
  AUTHOR: Ethan Brouwer + Codex draft
  LAST MODIFIED: 2026-02-12
*/

// ---------------------------------------------------------------------------
// 1) CONFIG
// ---------------------------------------------------------------------------
var cfg = {
  // Core assets
  mrdsAsset: "projects/metalminingpersonalcopy/assets/USGS_MRDS",
  gaul0Asset: "FAO/GAUL/2015/level0",
  demAsset: "USGS/SRTMGL1_003",

  // Optional fallback if geometry was not imported from CSV
  xField: "X",
  yField: "Y",
  nameField: "Name",
  descField: "description",

  // Analysis window
  startYear: 1984,
  endYear: ee.Number(ee.Date(Date.now()).get("year")).subtract(1), // last full year
  includeCurrentPartialYear: false,
  includeYearsWithNoImages: false,

  // Memory controls
  singleSiteId: null,     // e.g. "12345" for one-site debugging
  singleSiteName: null,   // exact match to Name field
  sitePartitionCount: 1,  // set >1 to split sites into partitions
  sitePartitionIndex: 0,  // 0-based partition index
  exportPerYear: true,    // strongly recommended for memory safety

  // Seasonal compositing window (El Salvador default: dry season Nov-Apr)
  seasonStartMonth: 11,
  seasonEndMonth: 4,

  // Quality controls
  cloudCoverMax: 70,
  minIlluminationCosine: 0.1,
  minValidPixelPct: 20,

  // Buffer distances
  buffersM: [1000, 2000],

  // NDVI classes aligned to Montero et al. (2024)
  ndviBareMax: 0.1,
  ndviSparseMax: 0.2,

  // Terrain correction
  applyTopoCorrection: false,
  topoSlopeMinDeg: 5,

  // Export config
  exportFolder: "GEE_exports",
  exportPrefixLong: "mrds_mine_disturbance_long",
  exportPrefixWide: "mrds_mine_disturbance_wide",

  // Runtime
  scale: 30
};

// Landsat C2 L2 collections used by analysis year.
var L5 = "LANDSAT/LT05/C02/T1_L2";
var L7 = "LANDSAT/LE07/C02/T1_L2";
var L8 = "LANDSAT/LC08/C02/T1_L2";
var L9 = "LANDSAT/LC09/C02/T1_L2";

var reflectanceMult = 0.0000275;
var reflectanceAdd = -0.2;
var targetBands = ["blue", "green", "red", "nir", "swir1", "swir2"];
var l57Bands = ["SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B7"];
var l89Bands = ["SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7"];

// ---------------------------------------------------------------------------
// 2) STUDY AREA + MRDS INGEST
// ---------------------------------------------------------------------------
var elsal = ee.FeatureCollection(cfg.gaul0Asset)
  .filter(ee.Filter.eq("ADM0_NAME", "El Salvador"));

Map.setOptions("SATELLITE");
Map.centerObject(elsal, 8);
Map.addLayer(elsal.style({color: "FFFFFF", width: 2, fillColor: "00000000"}), {}, "El Salvador");

function parseCoord(value) {
  var s = ee.String(value).trim();
  var hasNum = s.match("^-?\\d+(\\.\\d+)?$").size().gt(0);
  return ee.Algorithms.If(hasNum, ee.Number.parse(s), null);
}

function makePointFromXY(f) {
  var x = parseCoord(f.get(cfg.xField));
  var y = parseCoord(f.get(cfg.yField));
  var invalid = ee.Algorithms.If(
    ee.Algorithms.IsEqual(x, null),
    1,
    ee.Algorithms.If(ee.Algorithms.IsEqual(y, null), 1, 0)
  );
  var invalidFlag = ee.Number(invalid).eq(1);
  var valid = ee.Number(ee.Algorithms.If(invalidFlag, 0, 1));
  var safeX = ee.Number(ee.Algorithms.If(invalidFlag, 0, x));
  var safeY = ee.Number(ee.Algorithms.If(invalidFlag, 0, y));
  return ee.Feature(ee.Geometry.Point([safeX, safeY]), f.toDictionary())
    .set("valid_coord", valid);
}

var mrdsRaw = ee.FeatureCollection(cfg.mrdsAsset);

var mrds = mrdsRaw
  .filter(ee.Filter.notNull([cfg.xField, cfg.yField]))
  .filter(ee.Filter.neq(cfg.xField, ""))
  .filter(ee.Filter.neq(cfg.yField, ""))
  .map(makePointFromXY)
  .filter(ee.Filter.eq("valid_coord", 1))
  .filterBounds(elsal)
  .map(function (f) {
    var siteName = ee.String(ee.Algorithms.If(
      ee.Algorithms.IsEqual(f.get(cfg.nameField), null),
      "unknown_site",
      f.get(cfg.nameField)
    ));
    var siteId = ee.String(siteName).cat("_").cat(ee.String(f.id()));

    return f.set({
      site_id: siteId,
      site_name: siteName,
      prod_stage: "unknown",
      oper_type: "unknown",
      commodities: "unknown"
    });
  });

Map.addLayer(mrds.style({color: "FFAA00", pointSize: 5}), {}, "MRDS sites");
print("MRDS sites in El Salvador:", mrds.size());

function applySiteScope(fc) {
  if (cfg.singleSiteId !== null) {
    return fc.filter(ee.Filter.eq("site_id", cfg.singleSiteId));
  }
  if (cfg.singleSiteName !== null) {
    return fc.filter(ee.Filter.eq("site_name", cfg.singleSiteName));
  }
  if (cfg.sitePartitionCount > 1) {
    var withRand = fc.randomColumn("_part_rand", 1337);
    var width = 1 / cfg.sitePartitionCount;
    var lo = cfg.sitePartitionIndex * width;
    var hi = lo + width;
    return withRand
      .filter(ee.Filter.gte("_part_rand", lo))
      .filter(ee.Filter.lt("_part_rand", hi));
  }
  return fc;
}

mrds = applySiteScope(mrds);
print("MRDS sites after scope:", mrds.size());

function makeBuffers(fc, distanceM) {
  return fc.map(function (f) {
    return f.buffer(distanceM)
      .copyProperties(f)
      .set("buffer_m", distanceM);
  });
}

var buffers1k = makeBuffers(mrds, cfg.buffersM[0]);
var buffers2k = makeBuffers(mrds, cfg.buffersM[1]);
var siteBuffers = buffers1k.merge(buffers2k);

Map.addLayer(buffers1k.style({color: "00FFFF", width: 1, fillColor: "00000000"}), {}, "MRDS buffers 1 km", false);
Map.addLayer(buffers2k.style({color: "FF00FF", width: 1, fillColor: "00000000"}), {}, "MRDS buffers 2 km", false);

// ---------------------------------------------------------------------------
// 3) LANDSAT PREPROCESSING
// ---------------------------------------------------------------------------
function maskLandsatL2(img) {
  var qa = img.select("QA_PIXEL");
  var fill = qa.bitwiseAnd(1 << 0).neq(0);
  var dilated = qa.bitwiseAnd(1 << 1).neq(0);
  var cirrus = qa.bitwiseAnd(1 << 2).neq(0);
  var cloud = qa.bitwiseAnd(1 << 3).neq(0);
  var shadow = qa.bitwiseAnd(1 << 4).neq(0);
  var snow = qa.bitwiseAnd(1 << 5).neq(0);
  var sat = img.select("QA_RADSAT").eq(0);
  var clear = fill.or(dilated).or(cirrus).or(cloud).or(shadow).or(snow).eq(0);
  return img.updateMask(clear).updateMask(sat);
}

function getTerrainProducts(dem) {
  var terrain = ee.Terrain.products(dem);
  return {
    slopeRad: terrain.select("slope").multiply(Math.PI / 180),
    aspectRad: terrain.select("aspect").multiply(Math.PI / 180)
  };
}

// Simple illumination normalization. Keep toggle-enabled to compare with baseline.
function applyTopographicCorrection(img, terrain) {
  var sunAzDeg = ee.Number(ee.Algorithms.If(
    ee.Algorithms.IsEqual(img.get("SUN_AZIMUTH"), null),
    180,
    img.get("SUN_AZIMUTH")
  ));
  var sunElDeg = ee.Number(ee.Algorithms.If(
    ee.Algorithms.IsEqual(img.get("SUN_ELEVATION"), null),
    45,
    img.get("SUN_ELEVATION")
  ));
  var sunAz = sunAzDeg.multiply(Math.PI / 180);
  var sunEl = sunElDeg.multiply(Math.PI / 180);
  var sunZen = ee.Number(Math.PI / 2).subtract(sunEl);

  var cosZ = ee.Image.constant(sunZen.cos());
  var sinZ = ee.Image.constant(sunZen.sin());
  var cosS = terrain.slopeRad.cos();
  var sinS = terrain.slopeRad.sin();
  var cosAziDiff = ee.Image.constant(sunAz).subtract(terrain.aspectRad).cos();

  var ic = cosZ.multiply(cosS).add(sinZ.multiply(sinS).multiply(cosAziDiff))
    .max(cfg.minIlluminationCosine);
  var corrFactor = cosZ.divide(ic);
  var slopeMask = terrain.slopeRad.gte(cfg.topoSlopeMinDeg * Math.PI / 180);

  var corrected = img.select(["blue", "green", "red", "nir", "swir1", "swir2"])
    .multiply(corrFactor.where(slopeMask.not(), 1))
    .rename(["blue", "green", "red", "nir", "swir1", "swir2"]);

  return img.addBands(corrected, null, true).set("topo_corrected", 1);
}

function addNdvi(img) {
  var ndvi = img.normalizedDifference(["nir", "red"]).rename("NDVI");
  return img.addBands(ndvi);
}

function standardizeLandsatByBands(colId, sourceBands) {
  return ee.ImageCollection(colId)
    .filter(ee.Filter.lt("CLOUD_COVER", cfg.cloudCoverMax))
    .map(maskLandsatL2)
    .map(function (img) {
      var optical = img.select(sourceBands, targetBands)
        .multiply(reflectanceMult)
        .add(reflectanceAdd);
      return img.addBands(optical, null, true).select(targetBands);
    });
}

var dem = ee.Image(cfg.demAsset);
var terrain = getTerrainProducts(dem);

var landsatBase = standardizeLandsatByBands(L5, l57Bands)
  .merge(standardizeLandsatByBands(L7, l57Bands))
  .merge(standardizeLandsatByBands(L8, l89Bands))
  .merge(standardizeLandsatByBands(L9, l89Bands))
  .filterBounds(elsal);

var landsatPrepared = landsatBase
  .map(function (img) {
    return ee.Image(
      ee.Algorithms.If(cfg.applyTopoCorrection, applyTopographicCorrection(img, terrain), img.set("topo_corrected", 0))
    );
  })
  .map(addNdvi);

// ---------------------------------------------------------------------------
// 4) ANNUAL COMPOSITES + CLASSIFICATION
// ---------------------------------------------------------------------------
function analysisYears() {
  var endYear = ee.Number(
    ee.Algorithms.If(
      cfg.includeCurrentPartialYear,
      ee.Date(Date.now()).get("year"),
      cfg.endYear
    )
  );
  return ee.List.sequence(cfg.startYear, endYear);
}

function seasonWindow(year) {
  year = ee.Number(year);
  var crossYear = ee.Number(cfg.seasonStartMonth).gt(cfg.seasonEndMonth);

  var start = ee.Date.fromYMD(year, cfg.seasonStartMonth, 1);
  var endYear = ee.Number(ee.Algorithms.If(crossYear, year.add(1), year));
  var end = ee.Date.fromYMD(endYear, cfg.seasonEndMonth, 1).advance(1, "month");

  return {start: start, end: end};
}

function annualComposite(year) {
  year = ee.Number(year);
  var win = seasonWindow(year);
  var col = landsatPrepared.filterDate(win.start, win.end).filterBounds(elsal);
  var count = col.size();
  var empty = ee.Image.constant([0, 0, 0, 0, 0, 0, 0])
    .rename(["blue", "green", "red", "nir", "swir1", "swir2", "NDVI"])
    .updateMask(ee.Image(0));
  var comp = ee.Image(ee.Algorithms.If(count.gt(0), col.median(), empty))
    .set({
      year: year,
      start_date: win.start.format("YYYY-MM-dd"),
      end_date: win.end.format("YYYY-MM-dd"),
      image_count: count
    });
  return comp;
}

function classifyNdvi(ndvi) {
  var ndviClass = ee.Image(0)
    .where(ndvi.lt(cfg.ndviBareMax), 1)
    .where(ndvi.gte(cfg.ndviBareMax).and(ndvi.lt(cfg.ndviSparseMax)), 2)
    .where(ndvi.gte(cfg.ndviSparseMax), 3)
    .rename("ndvi_class");

  var miningSoil = ndviClass.eq(1).rename("mining_soil");
  var nonMiningSoil = ndviClass.neq(1).rename("non_mining_soil");

  return ndvi.addBands([ndviClass, miningSoil, nonMiningSoil]);
}

var years = analysisYears();
print("Configured year count:", years.size());

// ---------------------------------------------------------------------------
// 5) ZONAL STATS PER SITE x YEAR x BUFFER
// ---------------------------------------------------------------------------
function perFeatureYearStats(feature, annualImg) {
  var geom = feature.geometry();
  var ndvi = annualImg.select("NDVI");
  var ndviClass = annualImg.select("ndvi_class");
  var pixelHa = ee.Image.pixelArea().divide(10000);

  var ndviStats = ndvi.reduceRegion({
    reducer: ee.Reducer.mean()
      .combine(ee.Reducer.median(), null, true)
      .combine(ee.Reducer.stdDev(), null, true),
    geometry: geom,
    scale: cfg.scale,
    tileScale: 4,
    bestEffort: true,
    maxPixels: 1e10
  });

  var areaImg = ee.Image.cat([
    ndvi.mask().multiply(pixelHa).rename("valid_ha"),
    ndviClass.eq(1).multiply(pixelHa).rename("bare_ha"),
    ndviClass.eq(2).multiply(pixelHa).rename("sparse_ha"),
    ndviClass.eq(3).multiply(pixelHa).rename("veg_ha")
  ]);

  var areaStats = areaImg.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: geom,
    scale: cfg.scale,
    tileScale: 4,
    bestEffort: true,
    maxPixels: 1e10
  });

  var totalAreaHa = ee.Number(geom.area()).divide(10000);
  var validAreaHa = ee.Number(areaStats.get("valid_ha", 0));
  var bareHa = ee.Number(areaStats.get("bare_ha", 0));
  var sparseHa = ee.Number(areaStats.get("sparse_ha", 0));
  var vegHa = ee.Number(areaStats.get("veg_ha", 0));

  var validPct = ee.Algorithms.If(totalAreaHa.gt(0), validAreaHa.divide(totalAreaHa).multiply(100), 0);
  var barePct = ee.Algorithms.If(totalAreaHa.gt(0), bareHa.divide(totalAreaHa).multiply(100), 0);
  var miningSoilPct = barePct;
  var nonMiningSoilPct = ee.Algorithms.If(totalAreaHa.gt(0), sparseHa.add(vegHa).divide(totalAreaHa).multiply(100), 0);

  var qaFlag = ee.Algorithms.If(
    ee.Number(validPct).lt(cfg.minValidPixelPct),
    "low_valid_pixels",
    "ok"
  );

  var year = ee.Number(annualImg.get("year"));
  var siteId = ee.String(feature.get("site_id"));
  var bufferM = ee.Number(feature.get("buffer_m"));

  return ee.Feature(null, {
    site_id: siteId,
    site_name: feature.get("site_name"),
    prod_stage: feature.get("prod_stage"),
    commodities: feature.get("commodities"),
    year: year,
    buffer_m: bufferM,
    site_buffer_key: siteId.cat("_").cat(bufferM.format()),
    site_year_key: siteId.cat("_").cat(year.format()),
    mean_ndvi: ndviStats.get("NDVI_mean"),
    median_ndvi: ndviStats.get("NDVI_median"),
    ndvi_sd: ndviStats.get("NDVI_stdDev"),
    valid_px_pct: validPct,
    area_total_ha: totalAreaHa,
    area_bare_ha: bareHa,
    area_sparse_ha: sparseHa,
    area_veg_ha: vegHa,
    bare_pct: barePct,
    mining_soil_pct: miningSoilPct,
    non_mining_soil_pct: nonMiningSoilPct,
    topo_correction_applied: cfg.applyTopoCorrection,
    start_date: annualImg.get("start_date"),
    end_date: annualImg.get("end_date"),
    image_count: annualImg.get("image_count"),
    qa_flag: qaFlag
  });
}

function buildStatsForYear(year) {
  year = ee.Number(year);
  var comp = annualComposite(year);
  var classed = classifyNdvi(comp.select("NDVI"));
  var fullImg = comp.addBands(classed.select(["ndvi_class", "mining_soil", "non_mining_soil"]), null, true);
  var count = ee.Number(fullImg.get("image_count"));

  var perYear = siteBuffers.map(function (f) {
    return perFeatureYearStats(f, fullImg);
  });

  return ee.FeatureCollection(ee.Algorithms.If(
    cfg.includeYearsWithNoImages,
    perYear,
    ee.Algorithms.If(count.gt(0), perYear, ee.FeatureCollection([]))
  ));
}

function scopeTag() {
  function clean(s) {
    return String(s).replace(/[^A-Za-z0-9_]+/g, "_").slice(0, 40);
  }
  if (cfg.singleSiteName !== null) {
    return "site_" + clean(cfg.singleSiteName);
  }
  if (cfg.singleSiteId !== null) {
    return "siteid_" + clean(cfg.singleSiteId);
  }
  if (cfg.sitePartitionCount > 1) {
    return "part_" + clean(cfg.sitePartitionIndex) + "_of_" + clean(cfg.sitePartitionCount);
  }
  return "all_sites";
}

if (cfg.exportPerYear) {
  years.evaluate(function (yearList) {
    var tag = scopeTag();
    yearList.forEach(function (year) {
      var fc = buildStatsForYear(year);
      var y = String(year);
      var name = cfg.exportPrefixLong + "_" + tag + "_" + y;
      Export.table.toDrive({
        collection: fc,
        description: name,
        folder: cfg.exportFolder,
        fileNamePrefix: name,
        fileFormat: "CSV"
      });
    });
    if (yearList.length > 0) {
      var previewYear = yearList[0];
      var previewRow = buildStatsForYear(previewYear).first();
      print("CSV first-row preview (" + String(previewYear) + "):", previewRow);
    }
    print("Per-year export tasks queued:", yearList.length);
    print("Scope tag:", tag);
  });
}
