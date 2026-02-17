#!/usr/bin/env python3
"""
Generate dependency-free PNG trend charts from MRDS clean CSV output.

Designed for offline environments where matplotlib is unavailable.

Usage:
  python analysis/mrds_plot_png.py \
    --input-clean gee/groundwork/MRDSoutputs/Divisadero_Mine_2013_2025_clean.csv \
    --outdir analysis/figures/mrds_trends
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import struct
import zlib
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


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
    x: int
    y: float


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


def parse_int(value: str | None) -> int | None:
    f = parse_float(value)
    if f is None:
        return None
    return int(round(f))


def clean_slug(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s).strip("_")


def write_png_rgb(path: Path, width: int, height: int, rgb: bytearray) -> None:
    def chunk(tag: bytes, data: bytes) -> bytes:
        return (
            struct.pack(">I", len(data))
            + tag
            + data
            + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF)
        )

    raw = bytearray()
    row_bytes = width * 3
    for y in range(height):
        raw.append(0)  # filter type 0
        start = y * row_bytes
        raw.extend(rgb[start : start + row_bytes])

    png = bytearray()
    png.extend(b"\x89PNG\r\n\x1a\n")
    png.extend(
        chunk(
            b"IHDR",
            struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0),  # RGB
        )
    )
    png.extend(chunk(b"IDAT", zlib.compress(bytes(raw), level=9)))
    png.extend(chunk(b"IEND", b""))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(bytes(png))


class Canvas:
    def __init__(self, width: int, height: int, bg: tuple[int, int, int] = (255, 255, 255)) -> None:
        self.width = width
        self.height = height
        self.px = bytearray(width * height * 3)
        self.fill(bg)

    def fill(self, color: tuple[int, int, int]) -> None:
        r, g, b = color
        for i in range(0, len(self.px), 3):
            self.px[i] = r
            self.px[i + 1] = g
            self.px[i + 2] = b

    def set(self, x: int, y: int, color: tuple[int, int, int]) -> None:
        if x < 0 or y < 0 or x >= self.width or y >= self.height:
            return
        i = (y * self.width + x) * 3
        self.px[i] = color[0]
        self.px[i + 1] = color[1]
        self.px[i + 2] = color[2]

    def line(self, x0: int, y0: int, x1: int, y1: int, color: tuple[int, int, int]) -> None:
        dx = abs(x1 - x0)
        sx = 1 if x0 < x1 else -1
        dy = -abs(y1 - y0)
        sy = 1 if y0 < y1 else -1
        err = dx + dy
        while True:
            self.set(x0, y0, color)
            if x0 == x1 and y0 == y1:
                break
            e2 = 2 * err
            if e2 >= dy:
                err += dy
                x0 += sx
            if e2 <= dx:
                err += dx
                y0 += sy

    def rect_outline(
        self, x0: int, y0: int, x1: int, y1: int, color: tuple[int, int, int]
    ) -> None:
        self.line(x0, y0, x1, y0, color)
        self.line(x1, y0, x1, y1, color)
        self.line(x1, y1, x0, y1, color)
        self.line(x0, y1, x0, y0, color)

    def circle(self, cx: int, cy: int, r: int, color: tuple[int, int, int]) -> None:
        x = r
        y = 0
        err = 1 - x
        while x >= y:
            for px, py in (
                (x, y),
                (y, x),
                (-y, x),
                (-x, y),
                (-x, -y),
                (-y, -x),
                (y, -x),
                (x, -y),
            ):
                self.set(cx + px, cy + py, color)
            y += 1
            if err < 0:
                err += 2 * y + 1
            else:
                x -= 1
                err += 2 * (y - x) + 1


def ols_fit(points: list[Point]) -> tuple[float, float] | None:
    if len(points) < 2:
        return None
    xs = [p.x for p in points]
    ys = [p.y for p in points]
    xbar = sum(xs) / len(xs)
    ybar = sum(ys) / len(ys)
    den = sum((x - xbar) ** 2 for x in xs)
    if den == 0:
        return None
    slope = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys, strict=False)) / den
    intercept = ybar - slope * xbar
    return slope, intercept


def draw_series_png(
    output_png: Path,
    points: list[Point],
    width: int = 1100,
    height: int = 700,
) -> None:
    if len(points) < 2:
        return

    canvas = Canvas(width, height)
    margin_l = 90
    margin_r = 40
    margin_t = 40
    margin_b = 80
    plot_x0 = margin_l
    plot_y0 = margin_t
    plot_x1 = width - margin_r
    plot_y1 = height - margin_b
    plot_w = plot_x1 - plot_x0
    plot_h = plot_y1 - plot_y0

    years = [p.x for p in points]
    values = [p.y for p in points]
    min_x, max_x = min(years), max(years)
    min_y, max_y = min(values), max(values)
    if min_x == max_x or min_y == max_y:
        return

    y_pad = (max_y - min_y) * 0.08
    min_y -= y_pad
    max_y += y_pad

    # grid
    grid = (232, 232, 232)
    for i in range(1, 6):
        y = plot_y0 + int(i * plot_h / 6)
        canvas.line(plot_x0, y, plot_x1, y, grid)
    year_span = max_x - min_x
    x_ticks = min(year_span, 10)
    if x_ticks > 0:
        for i in range(1, x_ticks):
            x = plot_x0 + int(i * plot_w / x_ticks)
            canvas.line(x, plot_y0, x, plot_y1, grid)

    # axis
    axis = (60, 60, 60)
    canvas.rect_outline(plot_x0, plot_y0, plot_x1, plot_y1, axis)

    def xpix(x: int) -> int:
        return plot_x0 + int((x - min_x) / (max_x - min_x) * plot_w)

    def ypix(y: float) -> int:
        return plot_y1 - int((y - min_y) / (max_y - min_y) * plot_h)

    # series
    line_c = (36, 99, 235)
    pt_c = (18, 52, 120)
    pix = [(xpix(p.x), ypix(p.y)) for p in points]
    for i in range(1, len(pix)):
        canvas.line(pix[i - 1][0], pix[i - 1][1], pix[i][0], pix[i][1], line_c)
    for x, y in pix:
        canvas.circle(x, y, 2, pt_c)

    # trend line
    fit = ols_fit(points)
    if fit is not None:
        slope, intercept = fit
        y0 = slope * min_x + intercept
        y1 = slope * max_x + intercept
        canvas.line(xpix(min_x), ypix(y0), xpix(max_x), ypix(y1), (216, 27, 96))

    write_png_rgb(output_png, width, height, canvas.px)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate PNG trends from MRDS clean CSV.")
    parser.add_argument("--input-clean", required=True, help="Clean CSV from analysis/mrds_trends.py")
    parser.add_argument("--outdir", required=True, help="Output directory for PNGs")
    parser.add_argument(
        "--metrics",
        default=",".join(DEFAULT_METRICS),
        help="Comma-separated metric list to plot",
    )
    args = parser.parse_args()

    input_path = Path(args.input_clean)
    outdir = Path(args.outdir)
    metrics = [m.strip() for m in args.metrics.split(",") if m.strip()]

    with input_path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    groups: dict[tuple[str, str, str], list[Point]] = defaultdict(list)
    for row in rows:
        site_name = str(row.get("site_name", "")).strip()
        buffer_m = str(row.get("buffer_m", "")).strip()
        year = parse_int(row.get("year"))
        if site_name == "" or buffer_m == "" or year is None:
            continue
        for metric in metrics:
            val = parse_float(row.get(metric))
            if val is None:
                continue
            groups[(site_name, buffer_m, metric)].append(Point(year, val))

    manifest_rows: list[dict[str, str | int]] = []
    for (site_name, buffer_m, metric), points in sorted(groups.items()):
        points.sort(key=lambda p: p.x)
        if len(points) < 2:
            continue
        filename = (
            f"{clean_slug(site_name)}__{clean_slug(buffer_m)}m__{clean_slug(metric)}.png"
        )
        out_png = outdir / filename
        draw_series_png(out_png, points)
        manifest_rows.append(
            {
                "site_name": site_name,
                "buffer_m": buffer_m,
                "metric": metric,
                "n_points": len(points),
                "png_file": str(out_png),
            }
        )

    outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = outdir / "manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f, fieldnames=["site_name", "buffer_m", "metric", "n_points", "png_file"]
        )
        w.writeheader()
        for row in manifest_rows:
            w.writerow(row)

    print(f"Generated PNG files: {len(manifest_rows)}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
