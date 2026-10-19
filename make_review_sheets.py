#!/usr/bin/env python3
"""Place the individual vector PDFs at actual size on three labeled review sheets.

Optional dependency: reportlab. pypdf is included in the root requirements.
"""
from io import BytesIO
from pathlib import Path
import subprocess

from pypdf import PdfReader, PdfWriter, Transformation
from pypdf.generic import RectangleObject
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas
import plot_style as style

ROOT = Path(__file__).resolve().parent
OUT = ROOT/"validation/review"
OUT.mkdir(parents=True, exist_ok=True)
pdfmetrics.registerFont(TTFont("Arial", "/System/Library/Fonts/Supplemental/Arial.ttf"))


def review_sheet(name, size, panels, text):
    width, height = [v*72 for v in size]
    writer = PdfWriter()
    page = writer.add_blank_page(width=width, height=height)
    for stem, x, y in panels:
        if "/" in stem:
            source = ROOT/(stem+".pdf")
        elif stem.startswith("Fig4FG-"):
            source = ROOT/"BB-94-titration/results/combined"/(stem+".pdf")
        else:
            source = ROOT/"model-morphometry/results"/(stem+".pdf")
        panel = PdfReader(source).pages[0]
        if source.stem.startswith(("Fig4E-model-bud-count", "Fig4FG-BB94-titration")):
            # Show intentional off-canvas labels without changing the panel's
            # nominal dimensions, frame placement, or assembly coordinates.
            panel.cropbox = RectangleObject([-12, 0, float(panel.mediabox.right),
                                             float(panel.mediabox.top)])
        page.merge_transformed_page(panel, Transformation().translate(x*72, y*72))
    buffer = BytesIO()
    labels = canvas.Canvas(buffer, pagesize=(width, height))
    for caption, x, y, points in text:
        labels.setFont("Arial", points)
        labels.drawString(x*72, y*72, caption)
    labels.save()
    buffer.seek(0)
    page.merge_page(PdfReader(buffer).pages[0])
    path = OUT/(name+".pdf")
    with path.open("wb") as stream:
        writer.write(stream)
    style.outline_pdf_text(path)
    subprocess.run(["pdftoppm", "-r", "250", "-png", "-singlefile", str(path), str(OUT/name)], check=True)


review_sheet("Fig4A-5DE-acircularity-review", (5.6, 4.5), [
    ("Fig4A-acircularity-varying-thickness-tol-bright", .3, 2.6),
    ("Fig5D-acircularity-full-washout", .3, .65),
    ("Fig5E-acircularity-post-washout", 3.0, .65),
], [
    ("Model acircularity: panel mapping", .3, 4.15, 10),
    ("4A - Paul Tol Bright", .3, 3.72, 8),
    ("create_acircularity_plot_exp.py", .3, 2.43, 6),
    ("5D - Full washout time course", .3, 2.04, 8),
    ("5E - Post washout; reset at frame 75", 3.0, 2.04, 8),
    ("create_acircularity_plot_washout.py", .3, .42, 6),
    ("create_acircularity_plot_washoutcontrol.py", 3.0, .42, 6),
    ("Panels shown at their requested full sizes.", .3, .15, 6),
])

review_sheet("Fig4EFG-titration-review", (2.8, 4.05), [
    ("Fig4E-model-bud-count-tol-bright", .75, 2.55),
    ("Fig4FG-BB94-titration-tol-bright", .75, .65),
], [
    ("Figure 4E/F/G - Paul Tol Bright", .3, 3.78, 9),
    ("4E", .25, 3.10, 8), ("4F", .25, 1.85, 8), ("4G", .25, 1.22, 8),
    ("4E: 1.32 x 0.92 in; 4F/G: 1.32 x 1.51 in together.", .2, .43, 6),
    ("4E: 3 independent simulations per thickness.", .2, .29, 6),
    ("Circles: observations; black marks: mean +/- SD.", .2, .15, 6),
])
review_sheet("Fig4AB-review", (3.2, 3.85), [
    ("Fig4A-acircularity-varying-thickness-tol-bright", .75, 2.35),
    ("mesh-morphometry/results/Fig4B-asphericitychange_evolution_BB94_cleanmesh_plus-tol-bright", .75, .8),
], [
    ("Figure 4A / 4B - Paul Tol Bright", .3, 3.56, 10),
    ("4A - Model", .3, 3.40, 8),
    ("4B - Experimental", .3, 1.85, 8),
    ("Control / 1 / 5 uM match z = 0.001 / 0.005 / 0.01.", .3, .49, 6),
    ("4B: 0.5 pt individual curves at 30% opacity.", .3, .33, 6),
    ("Frame width: 40.80813 mm. Panels: 1.95 x 0.96 in.", .3, .17, 6),
])
print("Saved the labeled review sheets as PDF and PNG.")
