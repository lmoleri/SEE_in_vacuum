#!/usr/bin/env python3
"""
Digitize the 5 nm curve from Fig. 10 of the AIP Advances 095303 paper.

This script renders page 7, isolates the Fig. 10 plot region, and extracts
the topmost colored curve as the 5 nm fitted curve. It outputs a CSV with
energy (eV) and SEY values.
"""
import argparse
import os
import subprocess
import struct
import tempfile

import numpy as np
from PyPDF2 import PdfReader, PdfWriter


def _render_page(pdf_path: str, page_index: int, scale: int, tmp_dir: str) -> str:
    os.makedirs(tmp_dir, exist_ok=True)
    page_pdf = os.path.join(tmp_dir, f"fig10_page{page_index+1}.pdf")
    reader = PdfReader(open(pdf_path, "rb"))
    writer = PdfWriter()
    writer.add_page(reader.pages[page_index])
    with open(page_pdf, "wb") as f:
        writer.write(f)

    # Render to PNG via QuickLook
    subprocess.run(["qlmanage", "-t", "-s", str(scale), "-o", tmp_dir, page_pdf], check=True)
    base = os.path.basename(page_pdf) + ".png"
    # qlmanage writes to /private/tmp even if /tmp is requested
    for root in (tmp_dir, "/private/tmp"):
        candidate = os.path.join(root, base)
        if os.path.isfile(candidate):
            return candidate
    raise FileNotFoundError("Rendered PNG not found from qlmanage.")


def _png_to_bmp(png_path: str, bmp_path: str) -> str:
    subprocess.run(["sips", "-s", "format", "bmp", png_path, "--out", bmp_path], check=True)
    return bmp_path


def _load_bmp_rgb(bmp_path: str) -> np.ndarray:
    with open(bmp_path, "rb") as f:
        header = f.read(54)
        px_offset = struct.unpack("<I", header[10:14])[0]
        width = struct.unpack("<i", header[18:22])[0]
        height = struct.unpack("<i", header[22:26])[0]
        bpp = struct.unpack("<H", header[28:30])[0]
        if bpp != 32:
            raise ValueError(f"Expected 32bpp BMP, got {bpp}bpp")
        f.seek(px_offset)
        data = f.read()
    h = abs(height)
    w = width
    row_stride = w * 4
    arr = np.frombuffer(data, dtype=np.uint8).reshape((h, row_stride))
    # BGRA -> RGB
    b = arr[:, 0::4]
    g = arr[:, 1::4]
    r = arr[:, 2::4]
    rgb = np.stack([r, g, b], axis=-1)
    return rgb


def _find_plot_bbox(rgb: np.ndarray) -> tuple[int, int, int, int]:
    h, w, _ = rgb.shape
    r = rgb[:, :, 0].astype(np.int16)
    g = rgb[:, :, 1].astype(np.int16)
    b = rgb[:, :, 2].astype(np.int16)
    spread = np.maximum(np.maximum(r, g), b) - np.minimum(np.minimum(r, g), b)
    bright = (r + g + b) / 3.0
    colored = (spread > 40) & (bright > 50)
    colored[: h // 2, :] = False  # bottom half
    ys, xs = np.where(colored)
    if ys.size == 0:
        raise RuntimeError("Could not locate colored plot region.")
    pad_left = 220
    pad_right = 80
    pad_top = 80
    pad_bottom = 80
    x_min = max(xs.min() - pad_left, 0)
    x_max = min(xs.max() + pad_right, w - 1)
    y_min = max(ys.min() - pad_top, 0)
    y_max = min(ys.max() + pad_bottom, h - 1)
    return x_min, y_min, x_max, y_max


def _find_axes(rgb_crop: np.ndarray) -> tuple[int, int, int, int]:
    # Return (x_left, x_right, y_top, y_bottom) in crop coordinates
    r = rgb_crop[:, :, 0].astype(np.int16)
    g = rgb_crop[:, :, 1].astype(np.int16)
    b = rgb_crop[:, :, 2].astype(np.int16)
    bright = (r + g + b) / 3.0
    spread = np.maximum(np.maximum(r, g), b) - np.minimum(np.minimum(r, g), b)
    # Slightly looser threshold to include thin gray axis lines
    dark = (bright < 100) & (spread < 30)

    h, w = dark.shape
    # horizontal axis: search bottom 30%
    y_start = int(h * 0.68)
    row_counts = dark[y_start:, :].sum(axis=1)
    y_bottom_guess = y_start + int(np.argmax(row_counts))

    # vertical axis: pick the column with the strongest continuous dark run
    wlim = int(w * 0.7)
    col_counts = dark[:, :wlim].sum(axis=0)
    x_left = int(np.argmax(col_counts))

    # Use the longest contiguous dark segment in the axis column as the y-axis span
    col = dark[:, x_left]
    runs = []
    run_start = None
    for idx, val in enumerate(col):
        if val and run_start is None:
            run_start = idx
        elif not val and run_start is not None:
            runs.append((run_start, idx - 1))
            run_start = None
    if run_start is not None:
        runs.append((run_start, len(col) - 1))

    if runs:
        # Prefer the run that contains the bottom axis guess; otherwise the longest run.
        chosen = None
        for start, end in runs:
            if start <= y_bottom_guess <= end:
                chosen = (start, end)
                break
        if chosen is None:
            chosen = max(runs, key=lambda r: r[1] - r[0])
        y_top, y_bottom = int(chosen[0]), int(chosen[1])
    else:
        # Fallback to using the bottom guess and a reasonable top fraction
        y_bottom = y_bottom_guess
        y_top = int(h * 0.2)

    # Use the x-axis line (row near y_bottom) to find right bound
    row = dark[y_bottom, :]
    x_indices = np.where(row)[0]
    if x_indices.size:
        x_right = int(x_indices.max())
    else:
        x_right = int(w * 0.9)

    return x_left, x_right, y_top, y_bottom


def _detect_x_ticks(dark: np.ndarray, x_left: int, x_right: int, y_bottom: int) -> list[int]:
    # Measure tick lengths in a taller window above the x-axis (exclude the axis line)
    y0 = max(y_bottom - 60, 0)
    y1 = max(y_bottom - 2, 1)
    window = dark[y0:y1, x_left:x_right]
    if window.size == 0:
        return []
    col_counts = window.sum(axis=0)
    max_count = int(col_counts.max()) if col_counts.size else 0
    if max_count <= 0:
        return []
    thresh = max(6, int(max_count * 0.6))
    tick_cols = np.where(col_counts >= thresh)[0]
    if tick_cols.size == 0:
        return []
    ticks = []
    start = tick_cols[0]
    prev = tick_cols[0]
    for c in tick_cols[1:]:
        if c == prev + 1:
            prev = c
            continue
        ticks.append((start + prev) // 2)
        start = c
        prev = c
    ticks.append((start + prev) // 2)
    return [x_left + t for t in ticks]


def _detect_y_ticks(dark: np.ndarray, x_left: int, y_top: int, y_bottom: int) -> list[int]:
    # Measure tick lengths in a wider window to the right of the y-axis
    x0 = min(x_left + 2, dark.shape[1] - 1)
    x1 = min(x_left + 60, dark.shape[1])
    window = dark[y_top:y_bottom, x0:x1]
    if window.size == 0:
        return []
    row_counts = window.sum(axis=1)
    max_count = int(row_counts.max()) if row_counts.size else 0
    if max_count <= 0:
        return []
    thresh = max(6, int(max_count * 0.6))
    tick_rows = np.where(row_counts >= thresh)[0]
    if tick_rows.size == 0:
        return []
    ticks = []
    start = tick_rows[0]
    prev = tick_rows[0]
    for r in tick_rows[1:]:
        if r == prev + 1:
            prev = r
            continue
        ticks.append((start + prev) // 2)
        start = r
        prev = r
    ticks.append((start + prev) // 2)
    return [y_top + t for t in ticks]


def _rgb_to_hue(rgb: np.ndarray) -> np.ndarray:
    # rgb expected in uint8 [0,255]
    rgb_f = rgb.astype(np.float32) / 255.0
    r = rgb_f[:, :, 0]
    g = rgb_f[:, :, 1]
    b = rgb_f[:, :, 2]
    mx = np.maximum(np.maximum(r, g), b)
    mn = np.minimum(np.minimum(r, g), b)
    diff = mx - mn
    hue = np.zeros_like(mx)
    # Avoid division by zero
    mask = diff > 1e-6
    # Red is max
    idx = mask & (mx == r)
    hue[idx] = (60.0 * ((g[idx] - b[idx]) / diff[idx]) + 360.0) % 360.0
    # Green is max
    idx = mask & (mx == g)
    hue[idx] = 60.0 * ((b[idx] - r[idx]) / diff[idx]) + 120.0
    # Blue is max
    idx = mask & (mx == b)
    hue[idx] = 60.0 * ((r[idx] - g[idx]) / diff[idx]) + 240.0
    return hue


def _digitize_top_curve(
    rgb_crop: np.ndarray,
    axes,
    x_max_fraction: float,
    hue_range: tuple[float, float] | None = None,
    min_saturation: float = 0.25,
) -> tuple[np.ndarray, np.ndarray]:
    x_left, x_right, y_top, y_bottom = axes
    r = rgb_crop[:, :, 0].astype(np.int16)
    g = rgb_crop[:, :, 1].astype(np.int16)
    b = rgb_crop[:, :, 2].astype(np.int16)
    spread = np.maximum(np.maximum(r, g), b) - np.minimum(np.minimum(r, g), b)
    bright = (r + g + b) / 3.0
    colored = (spread > 40) & (bright > 50)
    if hue_range is not None:
        hue = _rgb_to_hue(rgb_crop)
        rgb_f = rgb_crop.astype(np.float32) / 255.0
        r_f = rgb_f[:, :, 0]
        g_f = rgb_f[:, :, 1]
        b_f = rgb_f[:, :, 2]
        mx = np.maximum(np.maximum(r_f, g_f), b_f)
        mn = np.minimum(np.minimum(r_f, g_f), b_f)
        sat = np.zeros_like(mx)
        nonzero = mx > 1e-6
        sat[nonzero] = (mx[nonzero] - mn[nonzero]) / mx[nonzero]
        hmin, hmax = hue_range
        if hmin <= hmax:
            hue_mask = (hue >= hmin) & (hue <= hmax)
        else:
            # wrap-around case (e.g. 350-20)
            hue_mask = (hue >= hmin) | (hue <= hmax)
        colored = colored & hue_mask & (sat >= min_saturation)

    x_max = x_left + int((x_right - x_left) * x_max_fraction)
    xs = []
    ys = []
    for x in range(x_left + 1, x_max):
        col = colored[y_top:y_bottom, x]
        if not np.any(col):
            continue
        y_idx = np.where(col)[0].min()  # topmost colored pixel
        xs.append(x)
        ys.append(y_top + y_idx)
    return np.array(xs), np.array(ys)


def main():
    parser = argparse.ArgumentParser(description="Digitize Fig. 10 5 nm curve from paper PDF.")
    parser.add_argument("--pdf", required=True, help="Path to 095303 paper PDF")
    parser.add_argument("--page", type=int, default=7, help="1-based page index with Fig. 10 (default: 7)")
    parser.add_argument("--scale", type=int, default=2400, help="QuickLook render size (default: 2400)")
    parser.add_argument("--x-min", type=float, default=0.0, help="X-axis minimum (eV)")
    parser.add_argument("--x-max", type=float, default=1500.0, help="X-axis maximum (eV)")
    parser.add_argument("--y-min", type=float, default=0.0, help="Y-axis minimum (SEY)")
    parser.add_argument("--y-max", type=float, default=4.0, help="Y-axis maximum (SEY)")
    parser.add_argument("--x-tick-step", type=float, default=200.0, help="X-axis tick step (eV)")
    parser.add_argument("--y-tick-step", type=float, default=0.5, help="Y-axis tick step (SEY)")
    parser.add_argument("--x-origin", type=float, default=0.0, help="X-axis origin value at leftmost tick")
    parser.add_argument("--y-origin", type=float, default=0.0, help="Y-axis origin value at bottommost tick")
    parser.add_argument("--tick-mapping", choices=["auto", "off"], default="auto",
                        help="Use detected tick marks for mapping (auto) or force axis range mapping (off).")
    parser.add_argument("--debug", action="store_true", help="Print debug info about detected ticks/axes.")
    parser.add_argument("--x-max-fraction", type=float, default=0.82,
                        help="Fraction of plot width to use (exclude right legend)")
    parser.add_argument("--hue-min", type=float, default=None,
                        help="Optional hue filter min (degrees, 0-360) to isolate a curve color.")
    parser.add_argument("--hue-max", type=float, default=None,
                        help="Optional hue filter max (degrees, 0-360) to isolate a curve color.")
    parser.add_argument("--hue-min-sat", type=float, default=0.25,
                        help="Minimum saturation for hue filtering (default: 0.25).")
    parser.add_argument("--output", required=True, help="Output CSV path")
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tmp_dir:
        png_path = _render_page(args.pdf, args.page - 1, args.scale, tmp_dir)
        bmp_path = os.path.join(tmp_dir, "page.bmp")
        _png_to_bmp(png_path, bmp_path)
        rgb = _load_bmp_rgb(bmp_path)

    x_min, y_min, x_max, y_max = _find_plot_bbox(rgb)
    rgb_crop = rgb[y_min:y_max + 1, x_min:x_max + 1, :]
    axes = _find_axes(rgb_crop)
    hue_range = None
    if args.hue_min is not None and args.hue_max is not None:
        hue_range = (args.hue_min, args.hue_max)
    xs, ys = _digitize_top_curve(
        rgb_crop,
        axes,
        args.x_max_fraction,
        hue_range=hue_range,
        min_saturation=args.hue_min_sat,
    )

    # Build dark mask for ticks detection
    r = rgb_crop[:, :, 0].astype(np.int16)
    g = rgb_crop[:, :, 1].astype(np.int16)
    b = rgb_crop[:, :, 2].astype(np.int16)
    bright = (r + g + b) / 3.0
    spread = np.maximum(np.maximum(r, g), b) - np.minimum(np.minimum(r, g), b)
    dark = (bright < 70) & (spread < 25)

    # Map to data coordinates using ticks if possible
    x_left, x_right, y_top, y_bottom = axes
    x_ticks = _detect_x_ticks(dark, x_left, x_right, y_bottom)
    y_ticks = _detect_y_ticks(dark, x_left, y_top, y_bottom)

    use_ticks = args.tick_mapping == "auto" and len(x_ticks) >= 4 and len(y_ticks) >= 4
    if args.debug:
        print(f"Detected x ticks: {len(x_ticks)}; y ticks: {len(y_ticks)}; use_ticks={use_ticks}")
    if use_ticks:
        x_ticks = sorted(x_ticks)
        y_ticks = sorted(y_ticks)
        x_diffs = np.diff(x_ticks)
        y_diffs = np.diff(y_ticks)
        x_step_pix = float(np.median(x_diffs)) if x_diffs.size else (x_right - x_left)
        y_step_pix = float(np.median(y_diffs)) if y_diffs.size else (y_bottom - y_top)
        x_scale = args.x_tick_step / x_step_pix if x_step_pix > 0 else (args.x_max - args.x_min) / (x_right - x_left)
        y_scale = args.y_tick_step / y_step_pix if y_step_pix > 0 else (args.y_max - args.y_min) / (y_bottom - y_top)
        x_origin_px = x_ticks[0]
        y_origin_px = y_ticks[-1]  # bottommost tick
        x_data = args.x_origin + (xs - x_origin_px) * x_scale
        y_data = args.y_origin + (y_origin_px - ys) * y_scale
    else:
        x_data = args.x_min + (xs - x_left) / (x_right - x_left) * (args.x_max - args.x_min)
        y_data = args.y_min + (y_bottom - ys) / (y_bottom - y_top) * (args.y_max - args.y_min)

    # Sort by x
    order = np.argsort(x_data)
    x_data = x_data[order]
    y_data = y_data[order]

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w") as f:
        f.write("energy_eV,sey\n")
        for x, y in zip(x_data, y_data):
            f.write(f"{x:.3f},{y:.6f}\n")

    print(f"Saved digitized curve: {args.output} ({len(x_data)} points)")


if __name__ == "__main__":
    main()
