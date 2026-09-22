"""Reusable, auditable Matplotlib export; does not choose statistical analyses.

Use configure_style(), figure_mm(), then save_figure_bundle(). Scientific plots
must be computed from recorded data. This module never fabricates observations,
confidence intervals, significance stars, or journal endorsement.
"""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import platform
import shutil
import xml.etree.ElementTree as ET

import matplotlib as mpl
from matplotlib import font_manager
import matplotlib.pyplot as plt
from matplotlib.text import Text
import numpy as np
import pandas as pd

# Okabe-Ito/Wong palette. Use shape, fill, or direct labels as well as colour.
COLORS = {
    "black": "#000000", "orange": "#E69F00", "sky_blue": "#56B4E9",
    "green": "#009E73", "yellow": "#F0E442", "blue": "#0072B2",
    "vermillion": "#D55E00", "purple": "#CC79A7", "neutral": "#767676",
}
WIDTH_MM = {"single": 89.0, "one_and_half": 120.0, "double": 183.0}


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def configure_style(font_size: float = 7.0) -> dict:
    """Set physical-size defaults; target-journal instructions override these.

    PDF embeds TrueType fonts. SVG preserves editable text and therefore depends
    on the named font being installed; it does NOT embed the font itself.
    """
    if not 5 <= font_size <= 7:
        raise ValueError("Nature research-style body text should be 5–7 pt.")
    available = {entry.name for entry in font_manager.fontManager.ttflist}
    font = next((name for name in ("Arial", "Helvetica", "DejaVu Sans") if name in available))
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": [font],
        "font.size": font_size, "axes.titlesize": font_size,
        "axes.labelsize": font_size, "xtick.labelsize": font_size,
        "ytick.labelsize": font_size, "legend.fontsize": font_size,
        "figure.titlesize": font_size, "text.color": "black",
        "axes.labelcolor": "black", "xtick.color": "black", "ytick.color": "black",
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.6, "axes.grid": False,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "xtick.direction": "out", "ytick.direction": "out",
        "lines.linewidth": 0.8, "lines.markersize": 4,
        "legend.frameon": False, "legend.handlelength": 1.1,
        "legend.borderaxespad": 0.0, "legend.labelspacing": 0.6,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white", "savefig.transparent": False,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "axes.unicode_minus": True, "savefig.dpi": 450,
    })
    return {"font_family": font, "font_size_pt": font_size,
            "font_file": font_manager.findfont(font),
            "preferred_font_available": font in ("Arial", "Helvetica")}


def figure_mm(width_mm: float = WIDTH_MM["double"], height_mm: float = 110.0, **kwargs):
    if width_mm <= 0 or height_mm <= 0:
        raise ValueError("Figure dimensions must be positive.")
    return plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4), **kwargs)


def text_layout_check(fig) -> dict:
    """Detect off-canvas text at final physical size; visual review still needed."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bbox = fig.bbox
    outside, small, large, overlaps = [], [], [], []
    visible = []
    for item in fig.findobj(match=Text):
        if not item.get_visible() or not item.get_text().strip():
            continue
        visible.append(item)
        b = item.get_window_extent(renderer)
        # 0.5 px tolerance for harmless renderer boundary rounding.
        if b.x0 < bbox.x0 - .5 or b.y0 < bbox.y0 - .5 or b.x1 > bbox.x1 + .5 or b.y1 > bbox.y1 + .5:
            outside.append(item.get_text())
        if item.get_fontsize() < 5:
            small.append(item.get_text())
        if item.get_fontsize() > 8:
            large.append(item.get_text())
    for i, first in enumerate(visible):
        a = first.get_window_extent(renderer)
        for second in visible[i + 1:]:
            b = second.get_window_extent(renderer)
            if a.overlaps(b):
                overlaps.append([first.get_text(), second.get_text()])
    result = {"visible_text_items": len(visible), "off_canvas_text": outside,
              "text_below_5_pt": small, "text_above_8_pt": large,
              "overlapping_text_boxes": overlaps,
              "manual_visual_review_required": True}
    if outside or small or large or overlaps:
        raise ValueError(f"Figure layout check failed: {result}")
    return result


def save_figure_bundle(fig, stem: Path, source_data: pd.DataFrame, caption: str,
                       inputs: list[Path], producer: Path, details: dict,
                       dpi: int = 450) -> dict:
    """Export data, caption, editable vectors, preview, and a provenance manifest.

    Source data must be the observations/statistics actually plotted. Caption
    and details must explain n, pairing, summaries, error bars and limitations.
    The exporter preserves requested physical dimensions (no bbox_inches=tight).
    """
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    if any(stem.parent.glob(stem.name + ".*")):
        raise FileExistsError("Refusing to overwrite an existing figure bundle.")
    if source_data.empty or not caption.strip():
        raise ValueError("Provide plotted source data and a scientific caption.")
    numeric = source_data.select_dtypes(include="number")
    if not np.isfinite(numeric.to_numpy()).all():
        raise ValueError("Non-finite plotted data; handle missingness explicitly.")
    checks = text_layout_check(fig)
    source_data.to_csv(stem.with_suffix(".source_data.csv"), index=False, lineterminator="\n")
    stem.with_suffix(".caption.md").write_text(caption.rstrip() + "\n", encoding="utf-8")
    shutil.copyfile(producer, stem.with_suffix(".producer.py"))
    shutil.copyfile(__file__, stem.with_suffix(".style.py"))
    fig.savefig(stem.with_suffix(".pdf"), metadata={"Title": details.get("title", stem.name)})
    fig.savefig(stem.with_suffix(".svg"), metadata={"Title": details.get("title", stem.name)})
    fig.savefig(stem.with_suffix(".png"), dpi=dpi)
    pdf = stem.with_suffix(".pdf").read_bytes()
    svg = ET.fromstring(stem.with_suffix(".svg").read_text(encoding="utf-8"))
    svg_texts = svg.findall(".//{http://www.w3.org/2000/svg}text")
    svg_images = svg.findall(".//{http://www.w3.org/2000/svg}image")
    checks.update({"pdf_has_embedded_truetype": b"/FontFile2" in pdf,
                   "pdf_has_type3_font": b"/Subtype /Type3" in pdf,
                   "svg_editable_text_count": len(svg_texts),
                   "svg_raster_image_count": len(svg_images)})
    if not checks["pdf_has_embedded_truetype"] or checks["pdf_has_type3_font"] or not svg_texts:
        raise ValueError("Vector/font export verification failed.")
    files = sorted(stem.parent.glob(stem.name + ".*"))
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "figure_dimensions_mm": (fig.get_size_inches() * 25.4).tolist(),
        "png_dpi": dpi, "plotted_rows": len(source_data),
        "python": platform.python_version(), "matplotlib": mpl.__version__,
        "numpy": np.__version__, "pandas": pd.__version__,
        "producer": {"path": str(producer.resolve()), "sha256": sha256(producer)},
        "style": {"path": str(Path(__file__).resolve()), "sha256": sha256(__file__)},
        "inputs": [{"path": str(p.resolve()), "sha256": sha256(p)} for p in inputs],
        "checks": checks, "details": details,
        "artifacts": [{"file": p.name, "sha256": sha256(p), "bytes": p.stat().st_size} for p in files],
    }
    stem.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest
