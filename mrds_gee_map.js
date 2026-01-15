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
var mrds = ee.FeatureCollection("users/your-username/USGS_MRDS_ElSalvador");

// CSV field names (edit to match your import)
var nameField = "name";
var commodityField = "commod1";

// ================= PARAMETERS =================
var bufferMeters = 200;

// ================= MAP SETUP =================
Map.centerObject(mrds, 9);
Map.setOptions("SATELLITE");

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

// ================= HELPERS =================
function safeGetString(f, field, fallback) {
  var val = f.get(field);
  return ee.Algorithms.If(val, ee.String(val), fallback);
}

function normalizeKey(s) {
  return ee.String(s).trim().toLowerCase();
}

// ================= POINT STYLING =================
var styledPoints = mrds.map(function (f) {
  var commodRaw = ee.String(safeGetString(f, commodityField, "unknown"));
  var commodKey = normalizeKey(commodRaw);
  var color = commodityColors.get(commodKey, "00FF00");
  return f.set("style", {
    color: color,
    pointSize: 6,
    width: 1,
    commod_norm: commodKey
  });
});

// ================= BUFFER CREATION =================
// Buffers preserve MRDS metadata for future classifier use
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

// Count features per normalized commodity
var counts = styledPoints.aggregate_histogram("commod_norm");
var keys = counts.keys().sort();

keys.evaluate(function (kList) {
  kList.forEach(function (k) {
    var color = commodityColors.get(k, "00FF00");
    var count = counts.get(k);
    count.evaluate(function (c) {
      var row = ui.Panel({
        layout: ui.Panel.Layout.Flow("horizontal"),
        style: {margin: "0 0 4px 0"}
      });
      row.add(ui.Label({
        style: {
          backgroundColor: "#" + color.getInfo(),
          padding: "6px",
          margin: "0 6px 0 0"
        }
      }));
      row.add(ui.Label(k + " (" + c + ")"));
      legend.add(row);
    });
  });
});

Map.add(legend);

// ================= USER NOTE =================
print(
  "Click a point or buffer and use the Inspector panel to view MRDS metadata " +
  "(commodity, ore type, deposit type, development status, etc.)."
);
