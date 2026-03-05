/*
  FILE: gee/groundwork/san_sebastian_upstream_downstream_water_simple.js
  PURPOSE: Minimal upstream/downstream water screening script for San Sebastian Mine.

  WHAT THIS SCRIPT DOES
  - Builds a recent Sentinel-2 composite.
  - Derives a conservative water mask from optical water indices plus JRC water.
  - Splits the stream corridor into one upstream control and three downstream bins.
  - Summarizes three water-only screening proxies:
      * TSS_PROXY   (red reflectance proxy)
      * NDTI        (discoloration / suspended sediment proxy)
      * RED_GREEN   (red/green ratio proxy)

  IMPORTANT LIMITATIONS
  - Screening only. This does not measure mercury or dissolved metals.
  - HydroSHEDS is used as coarse stream context, not precise channel geometry.
  - The downstream direction still depends on a user-set azimuth and should be
    verified visually in the GEE Code Editor before interpretation.
*/

// ================= USER SETTINGS =================
var siteName = "San Sebastian Mine (MRDS)";
var siteLon = -87.93002;
var siteLat = 13.6509;

var aoiRadiusMeters = 5000;
var recentMonths = 12;
var cloudMax = 60;
var useSeasonFilter = false;
var seasonMonths = [11, 12, 1, 2, 3, 4];

var downstreamAzimuthDeg = 110;
var streamAccThreshold = 50;
var streamBufferMeters = 90;
var outletSnapMeters = 1000;
var downstreamBinEdgesMeters = [0, 500, 1500, 3000];
var upstreamControlRangeMeters = [-1500, -500];

var scaleMeters = 30;
var minWaterPixelsPerZone = 5;
var noDataValue = -9999;

// ================= MAP SETUP =================
Map.setOptions("SATELLITE");

var sitePoint = ee.Geometry.Point([siteLon, siteLat]);
var aoi = sitePoint.buffer(aoiRadiusMeters);

Map.centerObject(sitePoint, 13);
Map.addLayer(sitePoint, {color: "FF0000"}, siteName, true);
Map.addLayer(aoi, {color: "FFFFFF"}, "AOI (5 km)", false);

// ================= DATE WINDOW =================
var recentEnd = ee.Date(Date.now());
var recentStart = recentEnd.advance(-recentMonths, "month");

print("Recent window", recentStart.format("YYYY-MM-dd"), recentEnd.format("YYYY-MM-dd"));
print("Downstream bins (m)", downstreamBinEdgesMeters);
print("Upstream control range (m)", upstreamControlRangeMeters);

// ================= SENTINEL-2 HELPERS =================
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

function scaleS2(img) {
  var optical = img.select(["B2", "B3", "B4", "B8", "B11"]).multiply(0.0001);
  return img.addBands(optical, null, true);
}

function addMonth(img) {
  return img.set("month", ee.Date(img.get("system:time_start")).get("month"));
}

function addWaterBands(img) {
  var green = img.select("B3");
  var red = img.select("B4");
  var nir = img.select("B8");
  var swir1 = img.select("B11");

  var ndwi = img.normalizedDifference(["B3", "B8"]).rename("NDWI");
  var mndwi = img.normalizedDifference(["B3", "B11"]).rename("MNDWI");
  var ndti = img.normalizedDifference(["B4", "B3"]).rename("NDTI");
  var redGreen = red.divide(green).rename("RED_GREEN");
  var tssProxy = red.multiply(1000).rename("TSS_PROXY");

  return img.addBands([ndwi, mndwi, ndti, redGreen, tssProxy]);
}

function buildS2Collection(start, end) {
  var col = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
    .filterBounds(aoi)
    .filterDate(start, end)
    .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", cloudMax))
    .map(maskS2)
    .map(scaleS2)
    .map(addMonth);

  if (useSeasonFilter) {
    col = col.filter(ee.Filter.inList("month", seasonMonths));
  }

  return col.map(addWaterBands);
}

// ================= WATER CONTEXT =================
var jrcWater = ee.Image("JRC/GSW1_4/GlobalSurfaceWater")
  .select("occurrence")
  .clip(aoi);
var permanentWater = jrcWater.gte(50).rename("PERMANENT_WATER");

function waterMask(img) {
  var opticalWater = img.select("NDWI").gt(0.12).or(img.select("MNDWI").gt(0.15));
  return opticalWater.and(permanentWater);
}

// ================= STREAM CONTEXT =================
var hydroAcc = ee.Image("WWF/HydroSHEDS/15ACC").clip(aoi);
var streamMask = hydroAcc.gte(streamAccThreshold).rename("STREAM_MASK");

var azimuthRad = downstreamAzimuthDeg * Math.PI / 180;
var ux = Math.sin(azimuthRad);
var uy = Math.cos(azimuthRad);

var lonLat = ee.Image.pixelLonLat();
var signedAlong = lonLat.select("longitude").subtract(siteLon).multiply(ux)
  .add(lonLat.select("latitude").subtract(siteLat).multiply(uy))
  .rename("SIGNED_ALONG");

var snapCircle = sitePoint.buffer(outletSnapMeters);
var snappedPointFc = streamMask.selfMask().reduceToVectors({
  geometry: snapCircle,
  scale: scaleMeters,
  geometryType: "centroid",
  maxPixels: 1e9
}).map(function(f) {
  return f.set("dist_m", f.geometry().distance(sitePoint, 1));
});

var snappedPointCount = snappedPointFc.size();
var snappedPoint = ee.Feature(snappedPointFc.sort("dist_m").first());
var snappedGeom = ee.Geometry(ee.Algorithms.If(
  snappedPointCount.gt(0),
  snappedPoint.geometry(),
  sitePoint
));

var sourceImage = ee.Image().toByte().paint(snappedGeom, 1).selfMask();
var streamCost = ee.Image.constant(1).updateMask(streamMask);
var streamDistance = streamCost.cumulativeCost(sourceImage, aoiRadiusMeters)
  .rename("STREAM_DISTANCE_M");

var downstreamCorridor = streamMask
  .and(signedAlong.gte(0))
  .and(streamDistance.lte(aoiRadiusMeters))
  .focal_max(streamBufferMeters, "circle", "meters")
  .selfMask()
  .rename("DOWNSTREAM_CORRIDOR");

var upstreamCorridor = streamMask
  .and(signedAlong.lt(0))
  .and(streamDistance.lte(Math.abs(upstreamControlRangeMeters[0])))
  .focal_max(streamBufferMeters, "circle", "meters")
  .selfMask()
  .rename("UPSTREAM_CORRIDOR");

Map.addLayer(streamMask.selfMask(), {palette: ["00FFFF"]}, "HydroSHEDS stream mask", false);
Map.addLayer(snappedGeom, {color: "FFFF00"}, "Snapped outlet point", true);
Map.addLayer(downstreamCorridor, {palette: ["FF6D00"]}, "Downstream corridor", true);
Map.addLayer(upstreamCorridor, {palette: ["2962FF"]}, "Upstream corridor", true);

// ================= COMPOSITE + QA =================
var recentCol = buildS2Collection(recentStart, recentEnd);
var recent = recentCol.median().clip(aoi);
var recentObsCount = recentCol.select("B4").count().clip(aoi).rename("RECENT_OBS_COUNT");
var currentWater = waterMask(recent).selfMask().rename("CURRENT_WATER");

var waterDownstream = currentWater.and(downstreamCorridor).selfMask().rename("WATER_DOWNSTREAM");
var waterUpstream = currentWater.and(upstreamCorridor).selfMask().rename("WATER_UPSTREAM");

print("Recent image count", recentCol.size());
print("Snapped stream candidates", snappedPointCount);

Map.addLayer(recent, {bands: ["B4", "B3", "B2"], min: 0.02, max: 0.30}, "Recent composite", false);
Map.addLayer(permanentWater.selfMask(), {palette: ["4FC3F7"]}, "Permanent water (JRC)", false);
Map.addLayer(currentWater, {palette: ["00BFFF"]}, "Current water mask", true);
Map.addLayer(waterDownstream, {palette: ["FF6D00"]}, "QA downstream water pixels", false);
Map.addLayer(waterUpstream, {palette: ["2962FF"]}, "QA upstream water pixels", false);
Map.addLayer(
  recentObsCount,
  {min: 0, max: 20, palette: ["2b2b2b", "f7f7f7", "00ff00"]},
  "QA observation count",
  false
);

Map.addLayer(
  recent.select("TSS_PROXY").updateMask(currentWater),
  {min: 20, max: 300, palette: ["08306b", "41b6c4", "ffffbf", "fdae61", "d73027"]},
  "Water TSS proxy",
  true
);
Map.addLayer(
  recent.select("NDTI").updateMask(currentWater),
  {min: -0.3, max: 0.4, palette: ["2166AC", "F7F7F7", "B2182B"]},
  "Water NDTI",
  false
);
Map.addLayer(
  recent.select("RED_GREEN").updateMask(currentWater),
  {min: 0.7, max: 2.0, palette: ["F7FBFF", "FDAE61", "D73027"]},
  "Water red/green ratio",
  false
);

// ================= ZONES =================
var downstreamBreaks = ee.List(downstreamBinEdgesMeters);
var emptyZoneGeom = ee.Geometry.Point([siteLon + 1, siteLat + 1]).buffer(1);

function maskToZoneFeature(maskImg, props) {
  var vectors = maskImg.selfMask().reduceToVectors({
    geometry: aoi,
    scale: scaleMeters,
    geometryType: "polygon",
    maxPixels: 1e9
  });
  var hasZone = vectors.size().gt(0);
  var zoneGeom = ee.Geometry(ee.Algorithms.If(hasZone, vectors.geometry(), emptyZoneGeom));

  return ee.Feature(zoneGeom, ee.Dictionary(props).set("has_zone", hasZone));
}

function zoneMaskFromFeature(feature) {
  feature = ee.Feature(feature);
  var zoneType = ee.String(feature.get("zone_type"));
  var inner = ee.Number(feature.get("inner_m"));
  var outer = ee.Number(feature.get("outer_m"));
  var hasZone = feature.get("has_zone");

  var zoneMask = ee.Image(ee.Algorithms.If(
    zoneType.compareTo("upstream").eq(0),
    streamDistance.gte(Math.abs(outer)).and(streamDistance.lt(Math.abs(inner))).and(upstreamCorridor),
    streamDistance.gte(inner).and(streamDistance.lt(outer)).and(downstreamCorridor)
  ));

  return ee.Image(ee.Algorithms.If(hasZone, zoneMask.selfMask(), ee.Image(0).selfMask()));
}

var downstreamZones = ee.FeatureCollection(
  ee.List.sequence(0, downstreamBreaks.length().subtract(2)).map(function(i) {
    i = ee.Number(i);
    var inner = ee.Number(downstreamBreaks.get(i));
    var outer = ee.Number(downstreamBreaks.get(i.add(1)));
    return maskToZoneFeature(
      streamDistance.gte(inner).and(streamDistance.lt(outer)).and(downstreamCorridor),
      {
        zone_id: ee.String("down_").cat(inner.format()).cat("_").cat(outer.format()),
        zone_type: "downstream",
        inner_m: inner,
        outer_m: outer,
        label: inner.format().cat("-").cat(outer.format()).cat(" m")
      }
    );
  })
);

var upstreamZone = maskToZoneFeature(
  streamDistance.gte(Math.abs(upstreamControlRangeMeters[1]))
    .and(streamDistance.lt(Math.abs(upstreamControlRangeMeters[0])))
    .and(upstreamCorridor),
  {
    zone_id: "upstream_control",
    zone_type: "upstream",
    inner_m: upstreamControlRangeMeters[0],
    outer_m: upstreamControlRangeMeters[1],
    label: "Upstream control"
  }
);

Map.addLayer(
  downstreamZones.style({color: "FFA500", fillColor: "00000000", width: 2}),
  {},
  "Downstream bins",
  true
);
Map.addLayer(
  ee.FeatureCollection([upstreamZone]).style({color: "2962FF", fillColor: "00000000", width: 2}),
  {},
  "Upstream control zone",
  true
);

print("Downstream bin availability", downstreamZones.aggregate_array("has_zone"));
print("Upstream control available", upstreamZone.get("has_zone"));

// ================= SUMMARIES =================
function safeNumber(value) {
  return ee.Number(ee.Algorithms.If(ee.Algorithms.IsEqual(value, null), 0, value));
}

function summarizeZone(feature, image, waterPixels, reachType) {
  var zoneOnlyMask = zoneMaskFromFeature(feature);
  var zoneWaterMask = ee.Image.constant(1)
    .updateMask(waterPixels)
    .updateMask(zoneOnlyMask);

  var waterCount = safeNumber(zoneWaterMask.reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: aoi,
    scale: scaleMeters,
    bestEffort: true,
    maxPixels: 1e9
  }).values().get(0));

  var stats = image.select(["TSS_PROXY", "NDTI", "RED_GREEN"])
    .updateMask(waterPixels)
    .updateMask(zoneOnlyMask)
    .reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: aoi,
      scale: scaleMeters,
      bestEffort: true,
      maxPixels: 1e9
    });

  var enoughWater = waterCount.gte(minWaterPixelsPerZone);

  return ee.Feature(feature.geometry(), feature.toDictionary())
    .set({
      reach_type: reachType,
      valid_pixel_count: waterCount,
      image_count: recentCol.size(),
      tss_proxy_mean: ee.Algorithms.If(enoughWater, stats.get("TSS_PROXY"), noDataValue),
      ndti_mean: ee.Algorithms.If(enoughWater, stats.get("NDTI"), noDataValue),
      red_green_mean: ee.Algorithms.If(enoughWater, stats.get("RED_GREEN"), noDataValue)
    });
}

var upstreamSummary = ee.FeatureCollection([
  summarizeZone(upstreamZone, recent, waterUpstream, "upstream")
]);

var downstreamSummary = downstreamZones.map(function(feature) {
  return summarizeZone(feature, recent, waterDownstream, "downstream");
});

var zoneSummaries = upstreamSummary.merge(downstreamSummary);
print("Zone summaries", zoneSummaries);

var tssChart = ui.Chart.feature.byFeature(zoneSummaries, "label", ["tss_proxy_mean"])
  .setChartType("ColumnChart")
  .setOptions({
    title: "Recent TSS proxy by reach",
    hAxis: {title: "Reach"},
    vAxis: {title: "Relative TSS proxy"},
    colors: ["#D84315"]
  });
print(tssChart);

var ndtiChart = ui.Chart.feature.byFeature(zoneSummaries, "label", ["ndti_mean"])
  .setChartType("ColumnChart")
  .setOptions({
    title: "Recent NDTI by reach",
    hAxis: {title: "Reach"},
    vAxis: {title: "NDTI"},
    colors: ["#8E24AA"]
  });
print(ndtiChart);

var redGreenChart = ui.Chart.feature.byFeature(zoneSummaries, "label", ["red_green_mean"])
  .setChartType("ColumnChart")
  .setOptions({
    title: "Recent red/green ratio by reach",
    hAxis: {title: "Reach"},
    vAxis: {title: "Red/green ratio"},
    colors: ["#2E7D32"]
  });
print(redGreenChart);

// ================= PANEL =================
var panel = ui.Panel({
  style: {position: "top-right", width: "330px", padding: "8px"}
});

panel.add(ui.Label({
  value: "San Sebastian Upstream/Downstream Water Screen",
  style: {fontWeight: "bold", fontSize: "13px"}
}));
panel.add(ui.Label("Map checks before interpretation:"));
panel.add(ui.Label("[x] Snapped outlet point is plausible"));
panel.add(ui.Label("[x] Upstream/downstream corridors follow the river"));
panel.add(ui.Label("[x] Water mask stays on the channel"));
panel.add(ui.Label("[x] Zone polygons align with visible water pixels"));
panel.add(ui.Label("Water markers in this script:"));
panel.add(ui.Label("[x] TSS proxy (red reflectance)"));
panel.add(ui.Label("[x] NDTI"));
panel.add(ui.Label("[x] Red/green ratio"));
panel.add(ui.Label({
  value: "Screening proxies only. Use field sampling to test pH, sulfate, conductivity, turbidity/TSS, and metals.",
  style: {fontSize: "11px", color: "B22222"}
}));

Map.add(panel);

print("Simple water screening script ready. Verify geometry and QA layers in GEE before exporting or making impact claims.");
