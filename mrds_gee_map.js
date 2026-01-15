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
// Import the MRDS KML manually in the GEE Code Editor
// and update the asset path below.
var mrds = ee.FeatureCollection("users/your-username/USGS_MRDS_ElSalvador");

// ================= PARAMETERS =================
var bufferMeters = 200;

// ================= MAP SETUP =================
Map.centerObject(mrds, 9);
Map.setOptions("SATELLITE");

// ================= COMMODITY COLOR MAP =================
var commodityColors = ee.Dictionary({
  "Gold": "FFD700",
  "Gold, Silver": "FFD700",
  "Silver, Gold": "FFD700",
  "Silver": "1E90FF",
  "Copper": "FF4500",
  "Lead": "8A2BE2",
  "Zinc": "00CED1",
  "Iron": "B22222",
  "Sulfur": "FFFF66",
  "Sulfur-Pyrite": "FFFF66",
  "Perlite": "A9A9A9"
});

// ================= POINT STYLING =================
var styledPoints = mrds.map(function (f) {
  var commod = ee.String(f.get("commod1"));
  var color = commodityColors.get(commod, "00FF00");
  return f.set("style", {
    color: color,
    pointSize: 6,
    width: 1
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

// ================= LABELS =================
var labels = mrds.map(function (f) {
  var name = ee.String(f.get("name", "Unnamed"));
  var commod = ee.String(f.get("commod1", "Unknown"));
  return f.set("label", name.cat(" | ").cat(commod));
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

Map.addLayer(
  labels.style({
    pointSize: 0,
    labelProperty: "label",
    fontSize: 10,
    textColor: "FFFFFF"
  }),
  {},
  "Site labels",
  false
);

// ================= USER NOTE =================
print(
  "Click a point or buffer and use the Inspector panel to view MRDS metadata " +
  "(commodity, ore type, deposit type, development status, etc.)."
);
