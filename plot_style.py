"""Shared manuscript typography, line weights and compact figure layout."""
from __future__ import annotations

import json
import io
from html import escape
from pathlib import Path
import shutil
import subprocess
import tempfile

import matplotlib as mpl
from matplotlib import font_manager

ROOT = Path(__file__).resolve().parent
FRAME_WIDTH = 0.6
DATA_WIDTH = 0.8
MEAN_WIDTH = 1.0
ERROR_WIDTH = 0.8
AXIS_SIZE = 7
TICK_SIZE = 6
LEGEND_SIZE = 6
ANNOTATION_SIZE = 7
AXIS_LABEL_PAD = 5
MESH_AXIS_LABEL_PAD = 2
TEX_POINT_FACTOR = 72.27 / 72
COMPACT_DOT_DIAMETER_SCALE = 1.2
COMPACT_CIRCLE_DIAMETER = 2.0
COMPACT_BUD_RATIO_CIRCLE_DIAMETER = 1.5
COMPACT_CIRCLE_EDGE_WIDTH = 0.35
COMPACT_POINT_ZORDER = 6
COMPACT_SUMMARY_ZORDER = 7


def compact_circle_style(collection, color):
    """Hollow circular observations, with summary marks drawn above them."""
    collection.set_facecolor("white")
    collection.set_edgecolor(color)
    collection.set_linewidth(COMPACT_CIRCLE_EDGE_WIDTH)
    collection.set_zorder(COMPACT_POINT_ZORDER)


def palettes():
    return json.loads((ROOT / "plot_palettes.json").read_text())


def titration_colors(name, doses):
    settings = palettes()
    if list(doses) != settings["titration_dose_order_um"]:
        raise ValueError("The selected titration palette requires the seven configured doses.")
    palette = settings["titration"][name]
    return [palette["control"], *palette["treatments"]]


def thickness_colors(name, z_values):
    """Look up model-thickness colors shared by panels 4A and 4B."""
    palette = palettes()["thickness"][name]
    mapping = dict(zip(palette["z_order"], palette["colors"], strict=True))
    try:
        return [mapping[z] for z in z_values]
    except KeyError as error:
        raise ValueError(f"No {name} color for model thickness {error.args[0]}") from error


def enlarge_compact_dots(collection):
    """Enlarge marker shapes while retaining the original swarm spacing."""
    from matplotlib.transforms import Affine2D
    transform = Affine2D().scale(COMPACT_DOT_DIAMETER_SCALE)
    collection.set_paths([path.transformed(transform) for path in collection.get_paths()])

AREV_PREAMBLE = (
    r"\usepackage{fontspec}\setmainfont{Arial}\setsansfont{Arial}"
    r"\usepackage{amsmath}\usepackage{arevmath}"
)


def rc(font="Arial", arev_math=False):
    font_manager.findfont(font_manager.FontProperties(family=font), fallback_to_default=False)
    result = {
        "font.family": font, "font.size": TICK_SIZE,
        "axes.labelsize": AXIS_SIZE, "axes.titlesize": AXIS_SIZE,
        "axes.labelpad": AXIS_LABEL_PAD,
        "xtick.labelsize": TICK_SIZE, "ytick.labelsize": TICK_SIZE,
        "legend.fontsize": LEGEND_SIZE, "legend.title_fontsize": LEGEND_SIZE,
        "axes.linewidth": FRAME_WIDTH, "xtick.major.width": FRAME_WIDTH,
        "ytick.major.width": FRAME_WIDTH, "xtick.minor.width": FRAME_WIDTH,
        "ytick.minor.width": FRAME_WIDTH, "lines.linewidth": DATA_WIDTH,
        "svg.fonttype": "path", "pdf.fonttype": 42, "text.usetex": False,
    }
    if arev_math:
        # XeLaTeX fontsize uses TeX points. Compensate so exported PDF text is
        # exactly 6 or 7 PostScript points, matching native Matplotlib exports.
        for key in ["font.size", "axes.labelsize", "axes.titlesize", "xtick.labelsize",
                    "ytick.labelsize", "legend.fontsize", "legend.title_fontsize"]:
            result[key] *= TEX_POINT_FACTOR
        result.update({"pgf.texsystem": "xelatex", "pgf.rcfonts": False,
                       "pgf.preamble": AREV_PREAMBLE,
                       "axes.labelpad": MESH_AXIS_LABEL_PAD})
    return result


def dimensions(key):
    """Full exported canvas in inches, including linked future panel sizes."""
    figures = json.loads((ROOT / "figure_sizes.json").read_text())["figures"]
    visited = set()
    while True:
        if key in visited:
            raise ValueError(f"Circular figure size link: {key}")
        visited.add(key)
        settings = figures[key]
        if "same_size_as" not in settings:
            break
        key = settings["same_size_as"]
    return [value / 25.4 for value in settings["full_mm"]]


def frame_width_inches(key):
    """Resolve an explicit or linked plotting-frame width, without a render."""
    figures = json.loads((ROOT / "figure_sizes.json").read_text())["figures"]
    visited = set()
    while True:
        if key in visited:
            raise ValueError(f"Circular frame width link: {key}")
        visited.add(key)
        settings = figures[key]
        if "frame_width_as" in settings:
            key = settings["frame_width_as"]
            continue
        width = (settings["frame_mm"][0] if settings.get("frame_mm")
                 else settings.get("frame_width_mm"))
        if width is None:
            raise ValueError(f"Reference frame width has not been specified: {key}")
        return width / 25.4


def frame_height_inches(key):
    """Resolve an explicit or linked plotting-frame height, without a render."""
    figures = json.loads((ROOT / "figure_sizes.json").read_text())["figures"]
    visited = set()
    while True:
        if key in visited:
            raise ValueError(f"Circular frame height link: {key}")
        visited.add(key)
        settings = figures[key]
        if "frame_height_as" in settings:
            key = settings["frame_height_as"]
            continue
        height = (settings["frame_mm"][1] if settings.get("frame_mm")
                  else settings.get("frame_height_mm"))
        if height is None:
            raise ValueError(f"Reference frame height has not been specified: {key}")
        return height / 25.4


def compact_layout(fig, ax, key=None, arev_math=False, size_inches=None,
                   vector_mathtext=False):
    """Keep the configured outside margin and fixed physical font/line sizes.

    Use an explicit full canvas, without bbox_inches='tight', so the exported
    dimensions equal the requested dimensions. Measure labels in the renderer
    that will export them. Exact frame dimensions can be supplied in the config.
    """
    config = json.loads((ROOT / "figure_sizes.json").read_text())
    settings = config["figures"][key] if key else {}
    if size_inches is not None:
        fig.set_size_inches(size_inches)
    elif key:
        fig.set_size_inches(dimensions(key))
    if arev_math:
        from matplotlib.backends.backend_pgf import FigureCanvasPgf
        FigureCanvasPgf(fig)
    elif vector_mathtext:
        # Raster text metrics round to pixels and can falsely reject a small
        # fixed frame. Measure native vector text in points for PDF/SVG export.
        fig.set_dpi(72)
    pad = settings.get("margin_points", config["margin_points"]) / 72
    for _ in range(6):
        if vector_mathtext:
            from matplotlib.backends.backend_svg import RendererSVG
            renderer = RendererSVG(*fig.bbox.size, io.StringIO())
        else:
            renderer = fig.canvas.get_renderer()
        frame = ax.get_window_extent(renderer)
        tight = ax.get_tightbbox(renderer)
        left = max(0, frame.x0 - tight.x0) / fig.dpi + pad
        right = max(0, tight.x1 - frame.x1) / fig.dpi + pad
        bottom = max(0, frame.y0 - tight.y0) / fig.dpi + pad
        top = max(0, tight.y1 - frame.y1) / fig.dpi + pad
        width, height = fig.get_size_inches()
        available_w, available_h = width - left - right, height - bottom - top
        if min(available_w, available_h) <= 0:
            raise ValueError(f"Figure {key} is too small for its fixed-size labels")
        frame_mm = settings.get("frame_mm")
        if frame_mm:
            frame_w, frame_h = [x / 25.4 for x in frame_mm]
            if frame_w > available_w + 1e-5 or frame_h > available_h + 1e-5:
                if vector_mathtext and _fit_fallback_labels(
                        ax, frame_w - available_w, frame_h - available_h):
                    continue
                raise ValueError(f"Requested frame for {key} does not fit its full figure")
            left += (available_w - frame_w) / 2
            bottom += (available_h - frame_h) / 2
            available_w, available_h = frame_w, frame_h
        elif "frame_height_mm" in settings or "frame_height_as" in settings:
            frame_h = frame_height_inches(key)
            if frame_h > available_h + 1e-5:
                if vector_mathtext and _fit_fallback_labels(ax, 0, frame_h - available_h):
                    continue
                raise ValueError(f"Requested frame height for {key} does not fit its full figure")
            bottom += (available_h - frame_h) / 2
            available_h = frame_h
        if not frame_mm and ("frame_width_mm" in settings or "frame_width_as" in settings):
            frame_w = frame_width_inches(key)
            if frame_w > available_w + 1e-5:
                if vector_mathtext and _fit_fallback_labels(ax, frame_w - available_w, 0):
                    continue
                raise ValueError(f"Requested frame width for {key} does not fit its full figure")
            left += (available_w - frame_w) / 2
            available_w = frame_w
        ax.set_position([left / width, bottom / height,
                         available_w / width, available_h / height])


def _fit_fallback_labels(ax, width_deficit, height_deficit):
    """Allow approximate font metrics to reduce label padding, never the frame."""
    changes = [(ax.yaxis, width_deficit), (ax.xaxis, height_deficit)]
    reductions = [(axis, max(0, deficit) * 72 + 0.02 if deficit > 0 else 0)
                  for axis, deficit in changes]
    if any(reduction > axis.labelpad for axis, reduction in reductions):
        return False
    for axis, reduction in reductions:
        axis.labelpad -= reduction
    return True


def compact_x_annotations(fig, ax, prefix, unit=None, arev_math=False):
    """Add 6 pt text beside tick anchors without participating in layout.

    Keep the first tilted tick's vertical anchor. Right-align the prefix and
    optional second-line unit to the left edge of the existing y tick labels.
    Text may extend beyond the canvas; no frame or tick is moved.
    """
    width, height = fig.get_size_inches() * 72
    first = ax.get_xticklabels()[0]
    _, y = first.get_transform().transform(first.get_position()) * 72 / fig.dpi
    if arev_math:
        renderer = fig.canvas.get_renderer()
    else:
        from matplotlib.backends.backend_svg import RendererSVG
        renderer = RendererSVG(*fig.bbox.size, io.StringIO())
    x = min(label.get_window_extent(renderer).x0 for label in ax.get_yticklabels()
            if label.get_visible() and label.get_text()) * 72 / fig.dpi
    fontsize = 6 * (TEX_POINT_FACTOR if arev_math else 1)
    artist = fig.text(x/width, y/height, prefix, ha="right", va="top",
                      family=mpl.rcParams["font.family"], fontsize=fontsize, in_layout=False, clip_on=False)
    artist.set_gid("compact-x-prefix")
    if unit:
        artist = fig.text(x/width, (y-7)/height, unit, ha="right", va="top",
                          family=mpl.rcParams["font.family"], fontsize=fontsize, in_layout=False, clip_on=False)
        artist.set_gid("compact-x-unit")


def compact_titration_layout(fig, ax, key, arev_math=False, size_inches=None):
    """Wider 6 pt titration panel, retaining ticks and omitting the x title."""
    fig.set_size_inches(size_inches or dimensions(key))
    if arev_math:
        from matplotlib.backends.backend_pgf import FigureCanvasPgf
        FigureCanvasPgf(fig)
    else:
        fig.set_dpi(72)
    width, height = fig.get_size_inches() * 72
    left, bottom, right, top = 22, 24, 3, 4
    if width <= left + right or height <= bottom + top:
        raise ValueError("The titration panel is too small for 6 pt labels.")
    ax.set_position([left / width, bottom / height,
                     (width-left-right) / width, (height-bottom-top) / height])
    fontsize = 6 * (TEX_POINT_FACTOR if arev_math else 1)
    ax.yaxis.label.set_fontsize(fontsize)
    ax.yaxis.label.set_ha("center")
    ax.yaxis.label.set_va("center")
    center_y = (bottom + (height-bottom-top)/2) / height if key == "fig4-model-buds" else 0.5
    ax.yaxis.set_label_coords(4 / width, center_y, transform=fig.transFigure)
    ax.set_xlabel("")
    if key == "fig4-model-buds":
        compact_x_annotations(fig, ax, r"$z$:", arev_math=arev_math)


def stacked_titration_layout(fig, axes):
    """Aligned F/G frames and a single row of concentration tick labels."""
    fig.set_size_inches(dimensions("fig4-titration-combined"))
    fig.set_dpi(72)
    width, height = fig.get_size_inches() * 72
    left, bottom, right, top, gap = 22, 22, 3, 4, 8
    frame_height = (height - bottom - top - gap) / 2
    for ax, y in zip(axes, [bottom + frame_height + gap, bottom]):
        ax.set_position([left/width, y/height, (width-left-right)/width, frame_height/height])
        ax.set_xlabel("")
        ax.yaxis.label.set_fontsize(6)
        ax.yaxis.label.set_ha("center")
        ax.yaxis.label.set_va("center")
        ax.yaxis.set_label_coords(4/width, (y+frame_height/2)/height, transform=fig.transFigure)
    axes[0].tick_params(axis="x", bottom=False, labelbottom=False)
    compact_x_annotations(fig, axes[1], "[BB-94]", unit="(µM)")


def allow_svg_overflow(stem):
    """Retain annotations outside the nominal canvas for figure assembly."""
    path = Path(stem).with_suffix(".svg")
    text = path.read_text()
    path.write_text(text.replace("<svg ", "<svg overflow='visible' ", 1))


def render_pdf_preview(stem, dpi=300):
    """Use the outlined PDF's exact layout for the PNG preview."""
    command = shutil.which("pdftoppm")
    if not command:
        raise RuntimeError("PNG previews require Poppler's pdftoppm.")
    subprocess.run([command, "-r", str(dpi), "-png", "-singlefile",
                    str(Path(stem).with_suffix(".pdf")), str(stem)], check=True)


def outline_pdf_text(path):
    """Convert embedded TeX/review fonts to vector paths without rasterizing."""
    command = shutil.which("gs")
    if not command:
        raise RuntimeError("Outlined TeX/review PDFs require Ghostscript's gs.")
    path = Path(path).resolve()
    with tempfile.TemporaryDirectory(prefix="tauber-outline-") as tmp:
        output = Path(tmp)/"outlined.pdf"
        result = subprocess.run([
            command, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE",
            "-sDEVICE=pdfwrite", "-dCompatibilityLevel=1.7",
            "-dNoOutputFonts", "-dAutoRotatePages=/None",
            "-sColorConversionStrategy=LeaveColorUnchanged",
            f"-sOutputFile={output}", str(path),
        ], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f"PDF outlining failed:\n{result.stdout[-6000:]}")
        path.write_bytes(output.read_bytes())


def save_native_figure(fig, stem, formats=("pdf", "svg", "png"), dpi=300,
                       transparent=False, allow_overflow=False):
    """Export native Matplotlib text as paths, retaining off-canvas labels.

    Normal's path effect draws the original glyph outlines through the vector
    renderer. Apply it only during export, after layout and swarm placement,
    and restore artists afterward. This also outlines native mathtext in PDF.
    """
    from matplotlib.patheffects import Normal
    from matplotlib.text import Text

    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    texts = fig.findobj(match=Text)
    effects = [(text, text.get_path_effects()) for text in texts]
    try:
        for text in texts:
            text.set_path_effects([Normal()])
        with mpl.rc_context({"svg.fonttype": "path"}):
            for extension in formats:
                if extension != "png" or "pdf" not in formats:
                    fig.savefig(stem.with_suffix("."+extension), dpi=dpi,
                                transparent=transparent)
    finally:
        for text, original in effects:
            text.set_path_effects(original)
    if "svg" in formats:
        path = stem.with_suffix(".svg")
        svg = path.read_text()
        for text in texts:
            if text.get_gid():
                group = f'<g id="{text.get_gid()}">'
                svg = svg.replace(group, group[:-1]+f' aria-label="{escape(text.get_text(), quote=True)}">')
        path.write_text(svg)
        if allow_overflow:
            allow_svg_overflow(stem)
    if "png" in formats and "pdf" in formats:
        render_pdf_preview(stem, dpi=dpi)


def save_arev_figure(fig, stem, dpi=300):
    """Export Arial text and original Arev math as vector outlines.

    XeLaTeX renders the original glyphs, then Ghostscript outlines the PDF.
    A separate dvisvgm PGF driver retains vector geometry and outlines SVG
    glyphs directly. PNG is rendered from the final outlined PDF.
    """
    needed = ["xelatex", "xdvipdfmx", "dvisvgm", "pdftoppm", "gs"]
    missing = [command for command in needed if not shutil.which(command)]
    if missing:
        raise RuntimeError("Missing font/export tools: " + ", ".join(missing)
                           + ". Install TeX Live, Poppler and Ghostscript, or use --text-rendering mathtext.")
    stem = Path(stem).resolve()
    stem.parent.mkdir(parents=True, exist_ok=True)

    def run(command, cwd):
        result = subprocess.run(command, cwd=cwd, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f"{command[0]} failed:\n{result.stdout[-6000:]}")

    with tempfile.TemporaryDirectory(prefix="tauber-figure-") as tmp:
        work = Path(tmp)
        fig.savefig(work / "figure.pgf", backend="pgf", transparent=True)
        for svg in [False, True]:
            name = "svg" if svg else "pdf"
            driver = r"\def\pgfsysdriver{pgfsys-dvisvgm.def}" if svg else ""
            (work / f"{name}.tex").write_text(
                r"\documentclass{standalone}" + driver + r"\usepackage{pgf}"
                + AREV_PREAMBLE + r"\begin{document}\input{figure.pgf}\end{document}")
            run(["xelatex", "-no-pdf", "-interaction=nonstopmode", "-halt-on-error",
                 f"{name}.tex"], work)
            if svg:
                run(["dvisvgm", "--no-fonts", "--bbox=papersize",
                     f"--output={stem.with_suffix('.svg')}", "svg.xdv"], work)
            else:
                run(["xdvipdfmx", "-q", "-d", "5", "-o", str(stem.with_suffix(".pdf")), "pdf.xdv"], work)
        if any(artist.get_gid() == "compact-x-prefix" for artist in fig.texts):
            allow_svg_overflow(stem)
        outline_pdf_text(stem.with_suffix(".pdf"))
        render_pdf_preview(stem, dpi=dpi)
