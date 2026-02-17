#!/usr/bin/env python3
"""
Generate dependency-free PNG trend charts from MRDS clean CSV output.

Designed for offline environments where matplotlib is unavailable.
Every chart includes:
- title
- x-axis label
- y-axis label
- numeric tick labels

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


FONT_5X7: dict[str, list[str]] = {
    " ": ["00000"] * 7,
    "-": ["00000", "00000", "00000", "11111", "00000", "00000", "00000"],
    "_": ["00000", "00000", "00000", "00000", "00000", "00000", "11111"],
    ".": ["00000", "00000", "00000", "00000", "00000", "00110", "00110"],
    ":": ["00000", "00110", "00110", "00000", "00110", "00110", "00000"],
    "(": ["00010", "00100", "01000", "01000", "01000", "00100", "00010"],
    ")": ["01000", "00100", "00010", "00010", "00010", "00100", "01000"],
    "%": ["11001", "11010", "00100", "01000", "10110", "00110", "00000"],
    "0": ["01110", "10001", "10011", "10101", "11001", "10001", "01110"],
    "1": ["00100", "01100", "00100", "00100", "00100", "00100", "01110"],
    "2": ["01110", "10001", "00001", "00010", "00100", "01000", "11111"],
    "3": ["11110", "00001", "00001", "01110", "00001", "00001", "11110"],
    "4": ["00010", "00110", "01010", "10010", "11111", "00010", "00010"],
    "5": ["11111", "10000", "11110", "00001", "00001", "10001", "01110"],
    "6": ["00110", "01000", "10000", "11110", "10001", "10001", "01110"],
    "7": ["11111", "00001", "00010", "00100", "01000", "01000", "01000"],
    "8": ["01110", "10001", "10001", "01110", "10001", "10001", "01110"],
    "9": ["01110", "10001", "10001", "01111", "00001", "00010", "11100"],
    "A": ["01110", "10001", "10001", "11111", "10001", "10001", "10001"],
    "B": ["11110", "10001", "10001", "11110", "10001", "10001", "11110"],
    "C": ["01110", "10001", "10000", "10000", "10000", "10001", "01110"],
    "D": ["11110", "10001", "10001", "10001", "10001", "10001", "11110"],
    "E": ["11111", "10000", "10000", "11110", "10000", "10000", "11111"],
    "F": ["11111", "10000", "10000", "11110", "10000", "10000", "10000"],
    "G": ["01110", "10001", "10000", "10111", "10001", "10001", "01110"],
    "H": ["10001", "10001", "10001", "11111", "10001", "10001", "10001"],
    "I": ["01110", "00100", "00100", "00100", "00100", "00100", "01110"],
    "J": ["00001", "00001", "00001", "00001", "10001", "10001", "01110"],
    "K": ["10001", "10010", "10100", "11000", "10100", "10010", "10001"],
    "L": ["10000", "10000", "10000", "10000", "10000", "10000", "11111"],
    "M": ["10001", "11011", "10101", "10001", "10001", "10001", "10001"],
    "N": ["10001", "11001", "10101", "10011", "10001", "10001", "10001"],
    "O": ["01110", "10001", "10001", "10001", "10001", "10001", "01110"],
    "P": ["11110", "10001", "10001", "11110", "10000", "10000", "10000"],
    "Q": ["01110", "10001", "10001", "10001", "10101", "10010", "01101"],
    "R": ["11110", "10001", "10001", "11110", "10100", "10010", "10001"],
    "S": ["01111", "10000", "10000", "01110", "00001", "00001", "11110"],
    "T": ["11111", "00100", "00100", "00100", "00100", "00100", "00100"],
    "U": ["10001", "10001", "10001", "10001", "10001", "10001", "01110"],
    "V": ["10001", "10001", "10001", "10001", "10001", "01010", "00100"],
    "W": ["10001", "10001", "10001", "10001", "10101", "11011", "10001"],
    "X": ["10001", "10001", "01010", "00100", "01010", "10001", "10001"],
    "Y": ["10001", "10001", "01010", "00100", "00100", "00100", "00100"],
    "Z": ["11111", "00001", "00010", "00100", "01000", "10000", "11111"],
}


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
        raw.append(0)
        start = y * row_bytes
        raw.extend(rgb[start : start + row_bytes])

    png = bytearray()
    png.extend(b"\x89PNG\r\n\x1a\n")
    png.extend(chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)))
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

    def rect_outline(self, x0: int, y0: int, x1: int, y1: int, color: tuple[int, int, int]) -> None:
        self.line(x0, y0, x1, y0, color)
        self.line(x1, y0, x1, y1, color)
        self.line(x1, y1, x0, y1, color)
        self.line(x0, y1, x0, y0, color)

    def circle(self, cx: int, cy: int, r: int, color: tuple[int, int, int]) -> None:
        x = r
        y = 0
        err = 1 - x
        while x >= y:
            for px, py in ((x, y), (y, x), (-y, x), (-x, y), (-x, -y), (-y, -x), (y, -x), (x, -y)):
                self.set(cx + px, cy + py, color)
            y += 1
            if err < 0:
                err += 2 * y + 1
            else:
                x -= 1
                err += 2 * (y - x) + 1

    def draw_char(self, x: int, y: int, ch: str, color: tuple[int, int, int], scale: int = 1) -> None:
        glyph = FONT_5X7.get(ch, FONT_5X7[" "])
        for gy, row in enumerate(glyph):
            for gx, bit in enumerate(row):
                if bit != "1":
                    continue
                for sy in range(scale):
                    for sx in range(scale):
                        self.set(x + gx * scale + sx, y + gy * scale + sy, color)

    def draw_text(self, x: int, y: int, text: str, color: tuple[int, int, int], scale: int = 1) -> None:
        cursor_x = x
        for ch in text:
            self.draw_char(cursor_x, y, ch, color, scale=scale)
            cursor_x += (5 * scale) + scale


def text_width(text: str, scale: int = 1) -> int:
    if not text:
        return 0
    return len(text) * (5 * scale + scale) - scale


def safe_label_text(text: str) -> str:
    t = text.upper().replace("|", " ").replace("/", " ")
    return "".join(ch if ch in FONT_5X7 else " " for ch in t)


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
    title: str,
    x_label: str,
    y_label: str,
    width: int = 1100,
    height: int = 700,
) -> None:
    if len(points) < 2:
        return

    canvas = Canvas(width, height)
    margin_l = 130
    margin_r = 40
    margin_t = 95
    margin_b = 120
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
    x_ticks = min(max(year_span, 2), 6)
    for i in range(1, int(x_ticks)):
        x = plot_x0 + int(i * plot_w / x_ticks)
        canvas.line(x, plot_y0, x, plot_y1, grid)

    # axis
    axis = (60, 60, 60)
    canvas.rect_outline(plot_x0, plot_y0, plot_x1, plot_y1, axis)

    def xpix(x: int) -> int:
        return plot_x0 + int((x - min_x) / (max_x - min_x) * plot_w)

    def ypix(y: float) -> int:
        return plot_y1 - int((y - min_y) / (max_y - min_y) * plot_h)

    # axis ticks + labels
    tick_c = (80, 80, 80)
    for i in range(int(x_ticks) + 1):
        t = i / int(x_ticks)
        xv = int(round(min_x + t * (max_x - min_x)))
        xp = xpix(xv)
        canvas.line(xp, plot_y1, xp, plot_y1 + 8, tick_c)
        lbl = safe_label_text(str(xv))
        canvas.draw_text(xp - text_width(lbl) // 2, plot_y1 + 16, lbl, tick_c, scale=1)

    y_ticks = 6
    for i in range(y_ticks + 1):
        t = i / y_ticks
        yv = max_y - t * (max_y - min_y)
        yp = plot_y0 + int(t * plot_h)
        canvas.line(plot_x0 - 8, yp, plot_x0, yp, tick_c)
        lbl = safe_label_text(f"{yv:.2f}")
        canvas.draw_text(plot_x0 - 10 - text_width(lbl), yp - 4, lbl, tick_c, scale=1)

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

    # labels
    txt_c = (35, 35, 35)
    title_txt = safe_label_text(title)
    x_txt = safe_label_text(x_label)
    y_txt = safe_label_text(y_label)
    canvas.draw_text((width - text_width(title_txt, 2)) // 2, 20, title_txt, txt_c, scale=2)
    canvas.draw_text((width - text_width(x_txt, 2)) // 2, height - 48, x_txt, txt_c, scale=2)
    canvas.draw_text(10, 16, y_txt, txt_c, scale=1)
    canvas.draw_text(10, 30, safe_label_text("TREND LINE = OLS"), (216, 27, 96), scale=1)

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
        filename = f"{clean_slug(site_name)}__{clean_slug(buffer_m)}m__{clean_slug(metric)}.png"
        out_png = outdir / filename
        metric_label = metric.replace("_", " ")
        draw_series_png(
            out_png,
            points,
            title=f"{site_name} | {metric_label} | {buffer_m} m buffer",
            x_label="Year",
            y_label=f"Value ({metric_label})",
        )
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
