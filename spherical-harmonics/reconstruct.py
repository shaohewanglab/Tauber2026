#!/usr/bin/env python3
"""Reconstruct Cartesian coordinate functions from a spherical surface map.

This is the SH integration/evaluation used for Figs. S3/S4, packaged without
local paths or the separate spherical-parameterization toolchain. Arrays use
x,y,z columns in micrometres; faces are zero-based triangles.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import sph_harm_y
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "example/241104-SMG5_f160_map.npz"


def polydata(vertices, faces):
    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(np.ascontiguousarray(vertices), deep=True))
    cells = vtk.vtkCellArray()
    packed = np.column_stack([np.full(len(faces), 3), faces]).astype(np.int64)
    cells.SetCells(len(faces), numpy_to_vtk(packed.ravel(), deep=True,
                                          array_type=vtk.VTK_ID_TYPE))
    mesh = vtk.vtkPolyData()
    mesh.SetPoints(points)
    mesh.SetPolys(cells)
    return mesh


def validate_mapping(vertices, faces, sphere_vertices):
    """Require a connected, closed genus-zero mesh and an unfolded sphere map."""
    if (vertices.ndim != 2 or vertices.shape[1] != 3
            or sphere_vertices.shape != vertices.shape
            or not np.isfinite(vertices).all()
            or not np.isfinite(sphere_vertices).all()):
        raise ValueError("vertices and sphere_vertices must be finite (N, 3) arrays")
    if (faces.ndim != 2 or faces.shape[1] != 3 or not len(faces)
            or not np.issubdtype(faces.dtype, np.integer)
            or faces.min() < 0 or faces.max() >= len(vertices)):
        raise ValueError("faces must contain valid zero-based integer triangles")
    radius_error = float(np.max(np.abs(np.linalg.norm(sphere_vertices, axis=1)-1)))
    if radius_error > 1e-8:
        raise ValueError("sphere_vertices must lie on the unit sphere")
    edges = np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]],
                                   faces[:, [2, 0]]]), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    graph = coo_matrix((np.ones(len(unique)), (unique[:, 0], unique[:, 1])),
                       shape=(len(vertices), len(vertices)))
    components = connected_components(graph, directed=False, return_labels=False)
    euler = len(vertices)-len(unique)+len(faces)
    if components != 1 or euler != 2 or not np.all(counts == 2):
        raise ValueError("Expected one closed manifold genus-zero surface, with no unused vertices")
    physical = vertices[faces]
    if np.any(np.linalg.norm(np.cross(physical[:, 1]-physical[:, 0],
                                     physical[:, 2]-physical[:, 0]), axis=1) == 0):
        raise ValueError("Physical mesh contains zero-area triangles")
    tri = sphere_vertices[faces]
    determinant = np.einsum("ij,ij->i", tri[:, 0], np.cross(tri[:, 1], tri[:, 2]))
    orientation = np.sign(np.median(determinant))
    if orientation == 0 or np.any(determinant*orientation <= 0):
        raise ValueError("Spherical map contains folded or degenerate triangles")
    denominator = (1 + np.einsum("ij,ij->i", tri[:, 0], tri[:, 1])
                     + np.einsum("ij,ij->i", tri[:, 1], tri[:, 2])
                     + np.einsum("ij,ij->i", tri[:, 2], tri[:, 0]))
    area = float((2*np.arctan2(determinant, denominator)).sum())
    if not np.isclose(abs(area), 4*np.pi, rtol=0, atol=1e-8):
        raise ValueError("Spherical triangles must cover the sphere once (area 4*pi)")
    return {"vertices": len(vertices), "faces": len(faces),
            "euler_characteristic": euler, "components": int(components),
            "folded_triangles": 0, "unit_radius_max_error": radius_error,
            "signed_spherical_area": area}


def load_mapping(path):
    with np.load(path, allow_pickle=False) as data:
        vertices = data["vertices"].astype(float)
        faces = data["faces"].copy()
        sphere = data["sphere_vertices"].astype(float)
    quality = validate_mapping(vertices, faces, sphere)
    return vertices, faces, sphere, quality


def sphere_grid(ntheta=128, nphi=256):
    u, weights = leggauss(ntheta)
    theta = np.arccos(u)
    phi = 2*np.pi*np.arange(nphi)/nphi
    sin_theta = np.sin(theta)[:, None]
    points = np.stack([sin_theta*np.cos(phi)[None, :],
                       sin_theta*np.sin(phi)[None, :],
                       np.broadcast_to(u[:, None], (ntheta, nphi))], axis=-1)
    return theta, weights, points


def inverse_samples(sphere, faces, vertices, points):
    """Intersect rays with the SPHERICAL mesh, then interpolate physical xyz."""
    locator = vtk.vtkStaticCellLocator()
    locator.SetDataSet(polydata(sphere, faces))
    locator.BuildLocator()
    t, sub_id, cell_id = vtk.mutable(0.), vtk.mutable(0), vtk.mutable(0)
    hit, barycentric = [0., 0., 0.], [0., 0., 0.]
    cell = vtk.vtkGenericCell()
    result = np.empty(points.shape)
    for i, direction in enumerate(points.reshape(-1, 3)):
        found = locator.IntersectWithLine((0., 0., 0.), direction*2, 1e-10,
                                         t, hit, barycentric, sub_id, cell_id, cell)
        if not found:
            raise ValueError(f"Inverse map missed spherical sample {i}")
        ids = faces[int(cell_id)]
        a, b = barycentric[:2]
        result.reshape(-1, 3)[i] = (vertices[ids[0]]*(1-a-b)
                                   + vertices[ids[1]]*a + vertices[ids[2]]*b)
    return result


def integrate(values, theta, weights, max_degree):
    """Compute a_lm = integral x(Omega) conjugate(Y_lm(Omega)) dOmega."""
    fourier = np.fft.fft(values, axis=1)*(2*np.pi/values.shape[1])
    coefficients = np.zeros((max_degree+1, max_degree+1, 3), complex)
    for ell in range(max_degree+1):
        for m in range(ell+1):
            polar = sph_harm_y(ell, m, theta, 0).real
            coefficients[ell, m] = np.sum(weights[:, None]*polar[:, None]
                                         * fourier[:, m], axis=0)
    return coefficients


def evaluate(coefficients, sphere, degrees):
    """Recover negative orders through conjugate symmetry for real coordinates."""
    theta = np.arccos(np.clip(sphere[:, 2], -1, 1))
    phi = np.arctan2(sphere[:, 1], sphere[:, 0])
    result, surfaces = np.zeros((len(sphere), 3)), {}
    for ell in range(max(degrees)+1):
        for m in range(ell+1):
            y = sph_harm_y(ell, m, theta, phi)
            # a_l,-m Y_l,-m + a_l,m Y_l,m = 2 Re(a_l,m Y_l,m).
            result += (1 if m == 0 else 2)*np.real(y[:, None]*coefficients[ell, m])
        if ell in degrees:
            surfaces[ell] = result.copy()
    return surfaces


def display_sphere(subdivisions=6):
    source = vtk.vtkPlatonicSolidSource()
    source.SetSolidTypeToOctahedron()
    source.Update()
    subdivide = vtk.vtkLinearSubdivisionFilter()
    subdivide.SetInputConnection(source.GetOutputPort())
    subdivide.SetNumberOfSubdivisions(subdivisions)
    subdivide.Update()
    mesh = subdivide.GetOutput()
    sphere = vtk_to_numpy(mesh.GetPoints().GetData()).astype(float)
    sphere /= np.linalg.norm(sphere, axis=1)[:, None]
    faces = vtk_to_numpy(mesh.GetPolys().GetData()).reshape(-1, 4)[:, 1:].copy()
    return sphere, faces


def outward_faces(vertices, faces):
    # A global winding reversal changes normals, never vertex geometry.
    tri = (vertices-vertices.mean(axis=0))[faces]
    volume6 = np.einsum("ij,ij->i", tri[:, 0], np.cross(tri[:, 1], tri[:, 2])).sum()
    return faces[:, [0, 2, 1]] if volume6 < 0 else faces.copy()


def write_ply(path, vertices, faces):
    """Standard xyz PLY, retaining double precision and zero-based indices."""
    with path.open("w", encoding="ascii", newline="\n") as stream:
        stream.write("ply\nformat ascii 1.0\ncomment coordinates x y z in micrometres\n")
        stream.write(f"element vertex {len(vertices)}\nproperty double x\n"
                     "property double y\nproperty double z\n")
        stream.write(f"element face {len(faces)}\n"
                     "property list uchar int vertex_indices\nend_header\n")
        np.savetxt(stream, vertices, fmt="%.17g")
        np.savetxt(stream, np.column_stack([np.full(len(faces), 3), faces]), fmt="%d")


def preview(vertices, faces, surfaces, display_faces, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    panels = [("Cleaned mesh", vertices, outward_faces(vertices, faces))]
    panels += [(f"n = {n}", v, outward_faces(v, display_faces)) for n, v in surfaces.items()]
    all_points = np.concatenate([v for _, v, _ in panels])
    center = (all_points.min(0)+all_points.max(0))/2
    span = np.ptp(all_points, axis=0).max()*1.12
    fig = plt.figure(figsize=(3*len(panels), 3.1), facecolor="white")
    light = np.array([-1., -1.5, 3.]); light /= np.linalg.norm(light)
    for i, (title, v, f) in enumerate(panels, 1):
        ax = fig.add_subplot(1, len(panels), i, projection="3d", proj_type="ortho")
        tri = (v-center)[f]
        normals = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
        normals /= np.maximum(np.linalg.norm(normals, axis=1)[:, None], 1e-30)
        brightness = .38+.58*np.maximum(normals@light, 0)
        colors = np.column_stack([brightness*.86, brightness*.91, brightness,
                                  np.ones(len(f))])
        ax.add_collection3d(Poly3DCollection(tri, facecolors=colors, edgecolors="none",
                                            linewidths=0, rasterized=True))
        ax.set(xlim=(-span/2, span/2), ylim=(-span/2, span/2), zlim=(-span/2, span/2))
        ax.set_box_aspect((1, 1, 1)); ax.view_init(elev=90, azim=-90, roll=180)
        ax.set_axis_off(); ax.set_title(title, fontsize=12, pad=-8)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=.93, wspace=0)
    fig.savefig(path, dpi=180, facecolor="white")
    plt.close(fig)


def analytic_check():
    theta, weights, grid = sphere_grid(64, 128)
    scale, center = np.array([60., 40., 20.]), np.array([10., 20., 30.])
    coefficients = integrate(grid*scale+center, theta, weights, 10)
    sphere, _ = display_sphere(3)
    errors = {n: float(np.max(np.abs(v-(sphere*scale+center))))
              for n, v in evaluate(coefficients, sphere, (1, 5, 10)).items()}
    if max(errors.values()) > 1e-9:
        raise AssertionError(f"Analytic ellipsoid reconstruction failed: {errors}")
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=ROOT/"results")
    parser.add_argument("--degrees", type=int, nargs="+", default=[5, 10, 20, 40])
    parser.add_argument("--ntheta", type=int, default=128)
    parser.add_argument("--nphi", type=int, default=256)
    parser.add_argument("--no-preview", action="store_true")
    parser.add_argument("--verify", action="store_true", help="Check bundled inputs and approved coefficients")
    args = parser.parse_args()
    degrees = sorted(set(args.degrees))
    if min(degrees) < 0 or args.ntheta <= max(degrees) or args.nphi <= 2*max(degrees):
        parser.error("Use nonnegative degrees, ntheta > max degree, and nphi > 2*max degree")
    if args.verify and (args.input.resolve() != DEFAULT_INPUT.resolve()
                        or args.ntheta != 128 or args.nphi != 256 or degrees != [5, 10, 20, 40]):
        parser.error("--verify uses the bundled input, 128 x 256 quadrature and degrees 5 10 20 40")
    if args.verify:
        for item in json.loads((ROOT/"source-data-manifest.json").read_text()):
            p = ROOT/item["path"]
            if hashlib.sha256(p.read_bytes()).hexdigest() != item["sha256"]:
                raise AssertionError(f"Example input changed: {p.name}")
    vertices, faces, sphere, quality = load_mapping(args.input)
    center = (vertices.min(0)+vertices.max(0))/2
    theta, weights, grid = sphere_grid(args.ntheta, args.nphi)
    print(f"Sampling inverse map at {args.ntheta} x {args.nphi} nodes...", flush=True)
    values = inverse_samples(sphere, faces, vertices-center, grid)
    coefficients = integrate(values, theta, weights, max(degrees))
    display, display_faces = display_sphere()
    surfaces = {n: v+center for n, v in evaluate(coefficients, display, degrees).items()}
    at_vertices = evaluate(coefficients, sphere, degrees)
    errors = {str(n): {"rms_um": float(np.sqrt(np.mean(np.sum((v-(vertices-center))**2, axis=1)))),
                       "max_um": float(np.max(np.linalg.norm(v-(vertices-center), axis=1)))}
              for n, v in at_vertices.items()}
    report = {"input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
              "coordinate_order": "x,y,z", "units": "micrometres",
              "quadrature": [args.ntheta, args.nphi], "degrees": degrees,
              "mapping_quality": quality, "center_xyz_um": center.tolist(),
              "errors_at_corresponding_mapping_vertices": errors,
              "error_definition": "Unweighted vertex correspondence error; not a closest-surface or Hausdorff distance"}
    if args.verify:
        with np.load(ROOT/"example/reference_coefficients.npz", allow_pickle=False) as ref:
            np.testing.assert_allclose(coefficients, ref["coefficients"], rtol=1e-11, atol=1e-10)
            np.testing.assert_array_equal(center, ref["center_xyz_um"])
            difference = float(np.max(np.abs(coefficients-ref["coefficients"])))
        reference = json.loads((ROOT/"example/provenance.json").read_text())
        for n, error in errors.items():
            for measure, actual in error.items():
                np.testing.assert_allclose(actual, reference["reference_errors"][n][measure],
                                           rtol=1e-11, atol=1e-10)
        analytic = analytic_check()
        report["verification"] = {"status": "passed", "coefficient_max_abs_difference_um": difference,
                                  "analytic_ellipsoid_max_abs_errors_um": analytic}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output_dir/"coefficients.npz", coefficients=coefficients,
                        center_xyz_um=center, coordinate_order="x,y,z", units="micrometres")
    for n, v in surfaces.items():
        f = outward_faces(v, display_faces)
        np.savez_compressed(args.output_dir/f"reconstruction_n{n:02d}.npz", vertices=v, faces=f)
        write_ply(args.output_dir/f"reconstruction_n{n:02d}.ply", v, f)
    if not args.no_preview:
        preview(vertices, faces, surfaces, display_faces, args.output_dir/"degree_comparison.png")
    (args.output_dir/"reconstruction_report.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"degrees": degrees, "errors_um": errors,
                      "verification": report.get("verification", {}).get("status", "not requested")}, indent=2))
    print(f"Saved meshes, coefficients and report to {args.output_dir}")


if __name__ == "__main__":
    main()
