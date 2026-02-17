#!/usr/bin/env python3
"""
Post-process GEE MRDS disturbance CSV exports into trend summaries.

Usage:
  python analysis/mrds_trends.py \
    --input gee/groundwork/MRDSoutputs/mrds_mine_disturbance_long_all_sites_1984_2025.csv
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import median


DEFAULT_METRICS = [
    "mean_ndvi",
    "mean_ndmi",
    "mean_ndbi",
    "mean_ndti",
    "mean_savi",
    "mean_bsi",
    "bare_pct",
    "mining_soil_pct",
    "non_mining_soil_pct",
    "valid_px_pct",
]


@dataclass
class Point:
    year: int
    value: float


def parse_float(value: str | None) -> float | None:
    if value is None:
        return None
    s = str(value).strip()
    if s == "":
        return None
    try:
        out = float(s)
    except ValueError:
        return None
    if math.isnan(out):
        return None
    return out


def parse_year(value: str | None) -> int | None:
    v = parse_float(value)
    if v is None:
        return None
    return int(round(v))


def ols_slope(points: list[Point]) -> float | None:
    if len(points) < 2:
        return None
    xs = [p.year for p in points]
    ys = [p.value for p in points]
    xbar = sum(xs) / len(xs)
    ybar = sum(ys) / len(ys)
    den = sum((x - xbar) ** 2 for x in xs)
    if den == 0:
        return None
    num = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys, strict=False))
    return num / den


def theil_sen_slope(points: list[Point]) -> float | None:
    if len(points) < 2:
        return None
    slopes: list[float] = []
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            dx = points[j].year - points[i].year
            if dx != 0:
                slopes.append((points[j].value - points[i].value) / dx)
    if not slopes:
        return None
    return float(median(slopes))


def trend_direction(slope: float | None, eps: float = 1e-4) -> str:
    if slope is None:
        return "insufficient_data"
    if slope > eps:
        return "increasing"
    if slope < -eps:
        return "decreasing"
    return "flat"


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute site/buffer trend summaries from GEE CSV.")
    parser.add_argument("--input", required=True, help="Input long-form CSV exported from GEE.")
    parser.add_argument(
        "--summary-out",
        default=None,
        help="Output summary CSV path. Default: <input>_trend_summary.csv",
    )
    parser.add_argument(
        "--clean-out",
        default=None,
        help="Output cleaned CSV path. Default: <input>_clean.csv",
    )
    parser.add_argument(
        "--min-years",
        type=int,
        default=8,
        help="Minimum data points required before reporting slopes (default: 8).",
    )
    parser.add_argument(
        "--metrics",
        default=",".join(DEFAULT_METRICS),
        help="Comma-separated metric columns to summarize.",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input CSV not found: {input_path}")

    base = input_path.with_suffix("")
    summary_out = Path(args.summary_out) if args.summary_out else Path(f"{base}_trend_summary.csv")
    clean_out = Path(args.clean_out) if args.clean_out else Path(f"{base}_clean.csv")
    metrics = [m.strip() for m in args.metrics.split(",") if m.strip()]

    with input_path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    if not rows:
        raise ValueError("Input CSV is empty.")

    clean_rows: list[dict[str, object]] = []
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)

    for row in rows:
        site_name = str(row.get("site_name", "")).strip()
        site_id = str(row.get("site_id", "")).strip()
        buffer_m = str(row.get("buffer_m", "")).strip()
        year = parse_year(row.get("year"))
        if year is None or site_name == "" or buffer_m == "":
            continue

        clean_row: dict[str, object] = {
            "site_name": site_name,
            "site_id": site_id,
            "buffer_m": buffer_m,
            "year": year,
            "image_count": parse_year(row.get("image_count")),
            "qa_flag": str(row.get("qa_flag", "")),
        }
        for m in metrics:
            clean_row[m] = parse_float(row.get(m))
        clean_rows.append(clean_row)
        grouped[(site_name, buffer_m)].append(clean_row)

    clean_rows.sort(key=lambda r: (str(r["site_name"]), str(r["buffer_m"]), int(r["year"])))

    summary_rows: list[dict[str, object]] = []
    for (site_name, buffer_m), group_rows in sorted(grouped.items()):
        by_metric: dict[str, list[Point]] = defaultdict(list)
        for r in group_rows:
            year = int(r["year"])
            for m in metrics:
                val = r.get(m)
                if isinstance(val, float):
                    by_metric[m].append(Point(year=year, value=val))

        for metric in metrics:
            pts = sorted(by_metric.get(metric, []), key=lambda p: p.year)
            n = len(pts)
            start_year = pts[0].year if pts else None
            end_year = pts[-1].year if pts else None
            start_val = pts[0].value if pts else None
            end_val = pts[-1].value if pts else None
            abs_change = (end_val - start_val) if (start_val is not None and end_val is not None) else None
            pct_change = None
            if start_val is not None and end_val is not None and start_val != 0:
                pct_change = ((end_val - start_val) / abs(start_val)) * 100

            slope_ols = ols_slope(pts) if n >= args.min_years else None
            slope_theil = theil_sen_slope(pts) if n >= args.min_years else None

            summary_rows.append(
                {
                    "site_name": site_name,
                    "buffer_m": buffer_m,
                    "metric": metric,
                    "n_years": n,
                    "start_year": start_year,
                    "end_year": end_year,
                    "start_value": start_val,
                    "end_value": end_val,
                    "abs_change": abs_change,
                    "pct_change": pct_change,
                    "ols_slope_per_year": slope_ols,
                    "theil_sen_slope_per_year": slope_theil,
                    "direction_ols": trend_direction(slope_ols),
                    "direction_theil_sen": trend_direction(slope_theil),
                }
            )

    write_csv(
        clean_out,
        clean_rows,
        ["site_name", "site_id", "buffer_m", "year", "image_count", "qa_flag", *metrics],
    )
    write_csv(
        summary_out,
        summary_rows,
        [
            "site_name",
            "buffer_m",
            "metric",
            "n_years",
            "start_year",
            "end_year",
            "start_value",
            "end_value",
            "abs_change",
            "pct_change",
            "ols_slope_per_year",
            "theil_sen_slope_per_year",
            "direction_ols",
            "direction_theil_sen",
        ],
    )

    print(f"Input rows: {len(rows)}")
    print(f"Clean rows written: {len(clean_rows)} -> {clean_out}")
    print(f"Trend summary rows written: {len(summary_rows)} -> {summary_out}")


if __name__ == "__main__":
    main()
