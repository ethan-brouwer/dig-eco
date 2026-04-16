/*
  FILE: gee/turbidity/rio_san_sebastian_2012_marn_points.js
  PURPOSE: Map 2012 MARN Rio San Sebastian water-quality stations with
  turbidity and cyanide values for comparison with remote-sensing turbidity
  analysis.

  SOURCE CSV
  /Users/EthanBrouwer/Documents/Projects/turb/rio_san_sebastian2012MARN.csv

  NOTES
  - Original coordinates were reported as DMS without hemisphere markers.
  - Coordinates are treated as latitude north and longitude west.
  - Earth Engine point order is [longitude, latitude].
  - Decimal longitudes are negative so the points plot in El Salvador.
  - Turbidez is reported as UNT in the source table.
  - Cianuro is reported as mg/L in the source table.
*/

var noDataValue = -9999;
var elSalvadorBounds = ee.Geometry.Rectangle([-90.2, 13.0, -87.6, 14.6], null, false);

// ================= STATION DATA =================
var stations = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point([-87.924722, 13.638056]), {
    station_id: "K1",
    latitude_dms: "13°38'17\"",
    longitude_dms: "87°55'29\"W",
    pH: 8.35,
    dissolved_oxygen_mg_l: 9.79,
    turbidity_unt: 4175,
    cyanide_mg_l: 0.002
  }),
  ee.Feature(ee.Geometry.Point([-87.923056, 13.638889]), {
    station_id: "K2",
    latitude_dms: "13°38'20\"",
    longitude_dms: "87°55'23\"W",
    pH: 8.4,
    dissolved_oxygen_mg_l: noDataValue,
    dissolved_oxygen_note: "Sin datos",
    turbidity_unt: 5565,
    cyanide_mg_l: 0.002
  }),
  ee.Feature(ee.Geometry.Point([-87.921389, 13.639444]), {
    station_id: "K3",
    latitude_dms: "13°38'22\"",
    longitude_dms: "87°55'17\"W",
    pH: 2.815,
    dissolved_oxygen_mg_l: 0.8,
    turbidity_unt: 7.4,
    cyanide_mg_l: 0.102
  }),
  ee.Feature(ee.Geometry.Point([-87.919444, 13.637778]), {
    station_id: "K4",
    latitude_dms: "13°38'16\"",
    longitude_dms: "87°55'10\"W",
    pH: 7.132,
    dissolved_oxygen_mg_l: 8.66,
    turbidity_unt: 13.3,
    cyanide_mg_l: 0.002
  }),
  ee.Feature(ee.Geometry.Point([-87.917222, 13.637222]), {
    station_id: "K5",
    latitude_dms: "13°38'14\"",
    longitude_dms: "87°55'02\"W",
    pH: 7.2815,
    dissolved_oxygen_mg_l: 9.86,
    turbidity_unt: 12.4,
    cyanide_mg_l: 0.002
  }),
  ee.Feature(ee.Geometry.Point([-87.915278, 13.636111]), {
    station_id: "K6",
    latitude_dms: "13°38'10\"",
    longitude_dms: "87°54'55\"W",
    pH: 7.4,
    dissolved_oxygen_mg_l: 8.22,
    turbidity_unt: 9.3,
    cyanide_mg_l: 0.002
  }),
  ee.Feature(ee.Geometry.Point([-87.915833, 13.634167]), {
    station_id: "K7",
    latitude_dms: "13°38'03\"",
    longitude_dms: "87°54'57\"W",
    pH: 7.515,
    dissolved_oxygen_mg_l: 8.39,
    turbidity_unt: 7.25,
    cyanide_mg_l: 0.002
  }),
  ee.Feature(ee.Geometry.Point([-87.914722, 13.6325]), {
    station_id: "K8",
    latitude_dms: "13°37'57\"",
    longitude_dms: "87°54'53\"W",
    pH: 7.69,
    dissolved_oxygen_mg_l: 8.04,
    turbidity_unt: 7.2,
    cyanide_mg_l: 0.002
  }),
  ee.Feature(ee.Geometry.Point([-87.913333, 13.631389]), {
    station_id: "K9",
    latitude_dms: "13°37'53\"",
    longitude_dms: "87°54'48\"W",
    pH: 7.791,
    dissolved_oxygen_mg_l: 8.26,
    turbidity_unt: 6.45,
    cyanide_mg_l: 0.002
  }),
  ee.Feature(ee.Geometry.Point([-87.913611, 13.629167]), {
    station_id: "K10",
    latitude_dms: "13°37'45\"",
    longitude_dms: "87°54'49\"W",
    pH: 7.818,
    dissolved_oxygen_mg_l: 8.38,
    turbidity_unt: 6.8,
    cyanide_mg_l: 0.002
  }),
  ee.Feature(ee.Geometry.Point([-87.912778, 13.627222]), {
    station_id: "K11",
    latitude_dms: "13°37'38\"",
    longitude_dms: "87°54'46\"W",
    pH: 8.052,
    dissolved_oxygen_mg_l: 9.6,
    turbidity_unt: 6.75,
    cyanide_mg_l: 0.002
  })
]);

var stationsWithFlags = stations.map(function(f) {
  var coords = f.geometry().coordinates();
  var turbidity = ee.Number(f.get("turbidity_unt"));
  var cyanide = ee.Number(f.get("cyanide_mg_l"));
  var dissolvedOxygen = ee.Number(f.get("dissolved_oxygen_mg_l"));
  return f.set({
    longitude_dd: coords.get(0),
    latitude_dd: coords.get(1),
    in_el_salvador_bbox: f.geometry().intersects(elSalvadorBounds, 1),
    dissolved_oxygen_chart: ee.Algorithms.If(
      dissolvedOxygen.eq(noDataValue),
      null,
      dissolvedOxygen
    ),
    turbidity_class: ee.Algorithms.If(
      turbidity.gte(1000),
      "extreme",
      ee.Algorithms.If(turbidity.gte(50), "high", "lower")
    ),
    cyanide_class: ee.Algorithms.If(cyanide.gt(0.002), "elevated", "reported_0.002")
  });
});

// ================= MAP =================
Map.setOptions("SATELLITE");
Map.centerObject(stationsWithFlags, 15);
Map.addLayer(elSalvadorBounds, {color: "FFFFFF"}, "El Salvador coordinate QA bounds", false);

var stationOutline = stationsWithFlags.style({
  color: "FFFFFF",
  pointSize: 10,
  pointShape: "circle",
  width: 2,
  fillColor: "00000000"
});

var turbidityStyled = stationsWithFlags.map(function(f) {
  var turbidity = ee.Number(f.get("turbidity_unt"));
  var color = ee.String(ee.Algorithms.If(
    turbidity.gte(1000),
    "7F0000",
    ee.Algorithms.If(turbidity.gte(50), "FC8D59", "2B83BA")
  ));
  return f.set("style", {
    color: "FFFFFF",
    fillColor: color,
    pointSize: 8,
    pointShape: "circle",
    width: 1
  });
}).style({styleProperty: "style"});

var cyanideStyled = stationsWithFlags.map(function(f) {
  var cyanide = ee.Number(f.get("cyanide_mg_l"));
  var color = ee.String(ee.Algorithms.If(cyanide.gt(0.002), "D73027", "1A9850"));
  return f.set("style", {
    color: "FFFFFF",
    fillColor: color,
    pointSize: 7,
    pointShape: "diamond",
    width: 1
  });
}).style({styleProperty: "style"});

Map.addLayer(stationOutline, {}, "MARN station outlines", true);
Map.addLayer(turbidityStyled, {}, "MARN turbidity classes", true);
Map.addLayer(cyanideStyled, {}, "MARN cyanide classes", false);

// ================= OUTPUTS =================
print("MARN Rio San Sebastian 2012 stations", stationsWithFlags);
print(
  "Stations outside El Salvador QA bounds",
  stationsWithFlags.filter(ee.Filter.eq("in_el_salvador_bbox", false))
);
print("Stations sorted by turbidity", stationsWithFlags.sort("turbidity_unt", false));
print("Stations sorted by cyanide", stationsWithFlags.sort("cyanide_mg_l", false));

var stationTable = ui.Chart.feature.byFeature(
  stationsWithFlags,
  "station_id",
  [
    "longitude_dd",
    "latitude_dd",
    "turbidity_unt",
    "cyanide_mg_l",
    "pH",
    "dissolved_oxygen_chart"
  ]
).setChartType("Table")
  .setOptions({
    allowHtml: true
  });
print(stationTable);

var turbidityChart = ui.Chart.feature.byFeature(
  stationsWithFlags,
  "station_id",
  ["turbidity_unt"]
).setChartType("ColumnChart")
  .setOptions({
    title: "Rio San Sebastian 2012 MARN Turbidity by Station",
    hAxis: {title: "Station"},
    vAxis: {title: "Turbidity (UNT)", logScale: true},
    legend: {position: "none"},
    colors: ["#7F0000"]
  });
print(turbidityChart);

var cyanideChart = ui.Chart.feature.byFeature(
  stationsWithFlags,
  "station_id",
  ["cyanide_mg_l"]
).setChartType("ColumnChart")
  .setOptions({
    title: "Rio San Sebastian 2012 MARN Cyanide by Station",
    hAxis: {title: "Station"},
    vAxis: {title: "Cyanide (mg/L)"},
    legend: {position: "none"},
    colors: ["#D73027"]
  });
print(cyanideChart);

var panel = ui.Panel({
  style: {
    width: "360px",
    position: "bottom-left",
    padding: "8px"
  }
});

panel.add(ui.Label({
  value: "Rio San Sebastian 2012 MARN",
  style: {fontWeight: "bold", fontSize: "16px"}
}));
panel.add(ui.Label("Turbidity circles: red = >=1000 UNT, orange = >=50 UNT, blue = <50 UNT."));
panel.add(ui.Label("Cyanide diamonds: red = >0.002 mg/L, green = reported 0.002 mg/L."));
panel.add(ui.Label("Highest turbidity: K2 = 5565 UNT; K1 = 4175 UNT."));
panel.add(ui.Label("Highest cyanide: K3 = 0.102 mg/L; all others reported 0.002 mg/L."));
panel.add(ui.Label("Click map points to inspect exact station values."));
Map.add(panel);

Export.table.toDrive({
  collection: stationsWithFlags,
  description: "rio_san_sebastian_2012_marn_water_quality_points",
  fileFormat: "CSV"
});
