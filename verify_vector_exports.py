#!/usr/bin/env python3
"""Verify that manuscript vector exports contain outlined text and no fonts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

from pypdf import PdfReader
from pypdf.generic import ContentStream

ROOT = Path(__file__).resolve().parent
TEXT_OPERATORS = {b"Tj", b"TJ", b"'", b'"'}
PATH_OPERATORS = {b"m", b"l", b"c", b"re"}


def audit_pdf(path):
    reader = PdfReader(path)
    counts = {"font_resources": 0, "text_operators": 0,
              "raster_images": 0, "vector_path_operators": 0}
    seen = set()

    def inspect(container):
        identity = id(container)
        if identity in seen:
            return
        seen.add(identity)
        resources = container.get("/Resources", {})
        if hasattr(resources, "get_object"):
            resources = resources.get_object()
        fonts = resources.get("/Font", {})
        if hasattr(fonts, "get_object"):
            fonts = fonts.get_object()
        counts["font_resources"] += len(fonts)
        stream = container.get_contents() if hasattr(container, "get_contents") else container
        if stream is not None:
            for _, operation in ContentStream(stream, reader).operations:
                counts["text_operators"] += operation in TEXT_OPERATORS
                counts["vector_path_operators"] += operation in PATH_OPERATORS
        objects = resources.get("/XObject", {})
        if hasattr(objects, "get_object"):
            objects = objects.get_object()
        for reference in objects.values():
            obj = reference.get_object()
            if obj.get("/Subtype") == "/Image":
                counts["raster_images"] += 1
            elif obj.get("/Subtype") == "/Form":
                inspect(obj)

    for page in reader.pages:
        inspect(page)
    assert counts["font_resources"] == 0, (path, counts)
    assert counts["text_operators"] == 0, (path, counts)
    assert counts["raster_images"] == 0, (path, counts)
    assert counts["vector_path_operators"] > 0, (path, counts)
    counts["page_sizes_points"] = [[float(page.mediabox.width), float(page.mediabox.height)]
                                  for page in reader.pages]
    return counts


def audit_svg(path):
    root = ET.parse(path).getroot()
    forbidden = {"text", "tspan", "textPath", "font", "font-face", "glyph",
                 "missing-glyph", "image", "foreignObject"}
    for node in root.iter():
        assert node.tag.split("}")[-1] not in forbidden, (path, node.tag)
    assert "@font-face" not in path.read_text(), path
    paths = sum(node.tag.endswith("}path") for node in root.iter())
    assert paths > 0, path
    return {"text_elements": 0, "font_definitions": 0, "raster_images": 0,
            "path_elements": paths,
            "size_points": [float(root.get(key).removesuffix("pt")) for key in ["width", "height"]]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-models", action="store_true")
    parser.add_argument("--include-review", action="store_true")
    args = parser.parse_args()
    folders = ["BB-94-titration/results/combined", "BB-94-washout/results/combined",
               "mesh-morphometry/results"]
    if not args.skip_models:
        folders.append("model-morphometry/results")
    checks = {}
    for folder in folders:
        for pdf in sorted((ROOT/folder).glob("*.pdf")):
            assert re.match(r"^Fig(?:S)?\d+[A-Z]+-", pdf.name), f"Missing plot ID: {pdf}"
            svg = pdf.with_suffix(".svg")
            assert svg.is_file(), svg
            record = {"pdf": audit_pdf(pdf), "svg": audit_svg(svg)}
            pdf_size = record["pdf"]["page_sizes_points"][0]
            svg_size = record["svg"]["size_points"]
            assert all(abs(a-b) < .001 for a, b in zip(pdf_size, svg_size)), (pdf, pdf_size, svg_size)
            checks[str(pdf.relative_to(ROOT))] = record
    expected = 6 if args.skip_models else 10
    assert len(checks) == expected, (len(checks), expected)
    reviews = {}
    if args.include_review:
        for pdf in sorted((ROOT/"validation/review").glob("*.pdf")):
            assert re.match(r"^Fig(?:S)?\d+[A-Z]+-", pdf.name), f"Missing plot IDs: {pdf}"
            reviews[str(pdf.relative_to(ROOT))] = audit_pdf(pdf)
        assert len(reviews) == 3, len(reviews)
    report = {"status": "passed", "plot_sets_verified": len(checks),
              "export_text": "Vector outlines; no installed fonts needed for viewing or placing",
              "font_resources": 0, "live_text_operators_or_elements": 0,
              "raster_images": 0, "plots": checks, "review_sheets": reviews}
    (ROOT/"validation/vector-outline-checks.json").write_text(json.dumps(report, indent=2)+"\n")
    print(f"Verified {len(checks)} outlined PDF/SVG pairs and {len(reviews)} review PDFs: no fonts, live text, or raster images.")


if __name__ == "__main__":
    main()
