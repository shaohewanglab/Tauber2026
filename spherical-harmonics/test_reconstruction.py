import unittest
import numpy as np
from scipy.special import sph_harm_y

import reconstruct as sh


class ReconstructionTests(unittest.TestCase):
    def test_offset_ellipsoid_is_recovered_at_degree_one(self):
        self.assertLess(max(sh.analytic_check().values()), 1e-9)

    def test_real_evaluation_matches_explicit_negative_orders(self):
        rng = np.random.default_rng(123)
        coefficients = np.zeros((6, 6, 3), complex)
        points = rng.normal(size=(25, 3))
        points /= np.linalg.norm(points, axis=1)[:, None]
        theta = np.arccos(points[:, 2]); phi = np.arctan2(points[:, 1], points[:, 0])
        expected = np.zeros((len(points), 3), complex)
        for ell in range(6):
            for m in range(ell+1):
                a = rng.normal(size=3)+1j*rng.normal(size=3)
                if m == 0:
                    a = a.real
                coefficients[ell, m] = a
                expected += sph_harm_y(ell, m, theta, phi)[:, None]*a
                if m:
                    expected += sph_harm_y(ell, -m, theta, phi)[:, None]*((-1)**m*a.conj())
        np.testing.assert_allclose(sh.evaluate(coefficients, points, [5])[5], expected.real,
                                   atol=1e-13)
        np.testing.assert_allclose(expected.imag, 0, atol=1e-13)

    def test_bundled_map_is_closed_and_folded_triangle_is_rejected(self):
        v, f, sphere, report = sh.load_mapping(sh.DEFAULT_INPUT)
        self.assertEqual(report["euler_characteristic"], 2)
        self.assertEqual(report["folded_triangles"], 0)
        broken = f.copy(); broken[0] = broken[0, [0, 2, 1]]
        with self.assertRaisesRegex(ValueError, "folded"):
            sh.validate_mapping(v, broken, sphere)


if __name__ == "__main__":
    unittest.main()
