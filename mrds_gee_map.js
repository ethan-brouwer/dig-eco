/*
  MRDS mine site visualization for El Salvador.

  PURPOSE
  - Display USGS MRDS mineral occurrence points
  - Style points by primary commodity
  - Add classifier-ready buffers (default 200 m)
  - Enable interactive inspection via GEE Inspector

  NOTE
  - MRDS locations represent reported occurrences or prospects,
    NOT confirmed active mining.
*/

// ================= USER INPUT =================
// Import the MRDS CSV/KML manually in the GEE Code Editor
// and update the asset path below.
var mrdsTable = ee.FeatureCollection("projects/metalminingpersonalcopy/assets/USGS_MRDS");

// CSV field names (edit to match your import)
var nameField = "Name";
var descField = "description";
var xField = "X";
var yField = "Y";

// ================= PARAMETERS =================
var bufferMeters = 200;

// ================= MAP SETUP =================
Map.setOptions("SATELLITE");

// ================= EL SALVADOR BORDER =================
var elsal = ee.FeatureCollection("FAO/GAUL/2015/level0")
  .filter(ee.Filter.eq("ADM0_NAME", "El Salvador"));
Map.addLayer(elsal.style({color: "FFFFFF", width: 2, fillColor: "00000000"}), {}, "El Salvador");
Map.centerObject(elsal, 8);

// ================= COMMODITY COLOR MAP =================
var commodityColors = ee.Dictionary({
  "gold": "FFD700",
  "gold, silver": "FFD700",
  "silver, gold": "FFD700",
  "silver": "1E90FF",
  "copper": "FF4500",
  "lead": "8A2BE2",
  "zinc": "00CED1",
  "iron": "B22222",
  "sulfur": "FFFF66",
  "sulfur-pyrite": "FFFF66",
  "perlite": "A9A9A9"
});

// Client-side mirror for legend colors
var commodityColorsClient = {
  "gold": "FFD700",
  "gold, silver": "FFD700",
  "silver, gold": "FFD700",
  "silver": "1E90FF",
  "copper": "FF4500",
  "lead": "8A2BE2",
  "zinc": "00CED1",
  "iron": "B22222",
  "sulfur": "FFFF66",
  "sulfur-pyrite": "FFFF66",
  "perlite": "A9A9A9"
};

// ================= HELPERS =================
function safeGetString(f, field, fallback) {
  var val = f.get(field);
  return ee.Algorithms.If(val, ee.String(val), fallback);
}

function normalizeKey(s) {
  return ee.String(s).trim().toLowerCase();
}

function extractCommodity(desc) {
  var s = ee.String(desc);
  var m = s.match("commod1</th><td>([^<]*)");
  return ee.Algorithms.If(m.size().gt(1), ee.String(m.get(1)), "unknown");
}

// ================= BUILD GEOMETRY =================
var mrds = mrdsTable.map(function (f) {
  var x = ee.Number(f.get(xField));
  var y = ee.Number(f.get(yField));
  var geom = ee.Geometry.Point([x, y]);
  return ee.Feature(geom, f.toDictionary());
});

// ================= POINT STYLING =================
var styledPoints = mrds.map(function (f) {
  var commodRaw = ee.String(extractCommodity(safeGetString(f, descField, "")));
  var commodKey = normalizeKey(commodRaw);
  var color = commodityColors.get(commodKey, "00FF00");
  return f.set("commod_norm", commodKey).set("style", {
    color: color,
    pointSize: 6,
    width: 1
  });
});

// ================= BUFFER CREATION =================
var buffers = mrds.map(function (f) {
  return f.buffer(bufferMeters)
    .copyProperties(f)
    .set({
      buffer_m: bufferMeters,
      source: "USGS_MRDS"
    });
});

// ================= ADD TO MAP =================
Map.addLayer(
  styledPoints.style({styleProperty: "style"}),
  {},
  "MRDS sites (by commodity)"
);

Map.addLayer(
  buffers.style({
    color: "FFFFFF",
    width: 2,
    fillColor: "00000000"
  }),
  {},
  "MRDS buffers (200 m)"
);

// ================= LEGEND =================
var legend = ui.Panel({
  style: {
    position: "top-left",
    padding: "8px 10px"
  }
});

legend.add(ui.Label({
  value: "MRDS commodities (count)",
  style: {fontWeight: "bold", fontSize: "12px"}
}));

var counts = styledPoints.aggregate_histogram("commod_norm");
counts.evaluate(function (countsObj) {
  var keys = Object.keys(countsObj || {}).sort();
  keys.forEach(function (k) {
    var color = commodityColorsClient[k] || "00FF00";
    var count = countsObj[k];
    var row = ui.Panel({
      layout: ui.Panel.Layout.Flow("horizontal"),
      style: {margin: "0 0 4px 0"}
    });
    row.add(ui.Label({
      style: {
        backgroundColor: "#" + color,
        padding: "6px",
        margin: "0 6px 0 0"
      }
    }));
    row.add(ui.Label(k + " (" + count + ")"));
    legend.add(row);
  });
});

Map.add(legend);

// ================= USER NOTE =================
print(
  "Click a point or buffer and use the Inspector panel to view MRDS metadata " +
  "(commodity, ore type, deposit type, development status, etc.)."
);
