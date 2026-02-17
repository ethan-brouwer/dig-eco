Phase I Package: MRDS Mine Environmental Characterization (El Salvador focus)

Purpose
- Share a single folder with the core Phase I assets for collaborators.
- This package is copy-based; original project files remain in their source locations.

Folder map
- data/raw/
  - Raw GEE export CSVs used for analysis.
- data/processed/
  - Cleaned trend tables and anomaly flags derived from raw exports.
- gee_scripts/
  - JavaScript scripts for Google Earth Engine workflows (including fast partitioned variant).
- outputs/
  - figures/: generated PNG trend plots.
  - gee_csv_examples/: earlier/export example CSVs from GEE groundwork outputs.
- supporting_info/
  - reports/: mine profile text reports (Divisadero, El Porvenir, Hormiguero).
  - tools/: Python scripts for trend summarization and PNG generation.

Primary run order for collaborators
1. Run GEE script from gee_scripts/ to export CSV.
2. Place export CSV in data/raw/.
3. Run supporting_info/tools/mrds_trends.py to create processed tables.
4. Run supporting_info/tools/mrds_plot_png.py to create figure outputs.
5. Review anomaly flags and mine reports for interpretation context.
