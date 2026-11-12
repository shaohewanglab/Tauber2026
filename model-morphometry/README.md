# Model acircularity and bud-count plots

From the repository root:

```bash
python model-morphometry/reproduce_model_plots.py
python model-morphometry/verify_model_plots.py
```

The root `python reproduce_plots.py` command includes these four model exports. PDF/SVG text is outlined; PNG previews are rendered from the PDFs. Rendering uses the root requirements, shared `plot_style.py`, Arial and TeX Live/Poppler/Ghostscript, as described in the root README.

| Panel | Input | Output stem | Full canvas (inches) |
| --- | --- | --- | --- |
| 4A | Final `ModelData/thickness_paper/*_exp_V00[1-3]_k12.npz` | `Fig4A-acircularity-varying-thickness-tol-bright` | 1.95 x 0.96 |
| 4E | `ModelData/bud_counts_summary.csv` | `Fig4E-model-bud-count-tol-bright` | 1.32 x 0.92 |
| 5D | Final three-run reference/washout archives | `Fig5D-acircularity-full-washout` | 1.60 x 1.24 |
| 5E | Final three-run reference/washout archives | `Fig5E-acircularity-post-washout` | 1.60 x 1.24 |

## Final Figure 4 data, 2026-11-10

The Figure 4 inputs were copied unchanged from the OneDrive paper folder `ModelDataV2/thickness_paper` and `ModelDataV2/bud_counts_summary_cusps.csv`. The CSV keeps its established repository name, `ModelData/bud_counts_summary.csv`. The six effective thicknesses are **0.001, 0.002, 0.003, 0.005, 0.008 and 0.01**. Each has three independent simulations, V001-V003. Only these 18 final coordinate archives are used for 4A; unselected z=0.015 and 0.02 final runs are not bundled.

4A shows the mean and sample SD (ddof=1) of all three runs. Each run's acircularity is baselined at its own frame 0 before summarizing. All runs contain 151 saved frames spaced by 0.01 in normalized model time, covering tau/tau_g=0-1.5. The integration step in filenames is not used as the saved-frame interval. The source contours have 400 nodes and 399 consecutive edges: area closes the contour with a chord; perimeter excludes that chord, preserving Justin's original calculation. `chi = perimeter/(2*sqrt(pi*area)) - 1`; the plotted value is `100*(chi-chi_baseline)`. `acircularity-curves.csv` contains every run, and `Fig4A-acircularity-summary.csv` contains the means, SDs and n.

4A/4E and the Illustrator showcase shapes use Paul Tol Bright purple, blue, cyan, green, yellow and red in ascending effective thickness order. The experimental 4B colors remain blue/green/red; their configured color associations now correspond to z=0.002/0.005/0.01. This is a shared visual association, not a new parameter fit. 4A retains the 4B frame width of 40.80813036 mm, y=0-110%, a two-column legend, and right-aligned endpoint tick 1.5.

4E contains all 18 supplied bud counts at tau/tau_g=1.5, with hollow 2 pt circles (0.35 pt colored outlines), black mean/sample SD and five two-sided single-step Dunnett comparisons versus z=0.001. The normal/common-variance approximation and small sample size remain stated in the statistics reports; no counts or variances are imputed. An unbudded shape counts as one. Adjusted P values in ascending comparison-z order are 0.3356949191, 0.0846301915, 0.00006922114055, 0.00000004304914446 and 0.00000004304914446; annotations are ns, ns, ****, ****, ****. No additional Holm correction is applied to this family.

The independent-run Spearman trend is recalculated from all 18 counts (rho=-0.9451620167, raw permutation P=0.000001). Root `write_figure4_trends.py` applies the existing separate Holm correction across E/F/G trends (adjusted P=0.000003). Comparisons, counts, mean/SD, assumptions, source hashes and trends are in `results/`. See [internal data history](../INTERNAL_NOTES.md).

## Final Figure 5 data, 2026-11-12

5D/E use the three supplied `ModelDataV2/washout/*_exp_V00[1-3]_washout.npz`
archives and matching final constant-thickness runs at z=0.001 and z=0.015.
The three z=0.001 references are shared with 4A; the three z=0.015 references
are included for 5D. The six newly bundled archives are copied without changing
any bytes, and both source manifests record their supplied paths and hashes.

Each group shows the mean and sample SD (ddof=1) of three independent runs.
Each run is baselined separately: frame 0 for 5D and frame 75 for 5E. Saved
frames span 0-150 at intervals of 0.01 normalized time. 5E therefore spans
0-0.75 after resetting time and acircularity at tau/tau_g=0.75. Calculations
retain the chord-closed area and open-edge perimeter used in 4A. Source
integration steps in filenames are not the saved-frame time interval.

`acircularity-curves.csv` records every run. `Fig5DE-acircularity-summary.csv`
records the five group mean/SD trajectories and n=3 at every displayed frame.
The verifier independently recomputes the node areas, edge lengths, per-run
baseline changes, all five means and sample SD bands. Plot colors, sizes,
labels, limits and line widths retain the established Figure 5 presentation.

## Legacy inputs and verification

The older eight coordinate archives remain bundled only for the frozen historical `validation/original-curves.json` fixture. The four original scripts are unchanged in `original_scripts/`. The October bud-count table is preserved under `validation/history`. Both source manifests record current input hashes and their original project paths.

Historical reference verification reruns the supplied scripts against the legacy inputs and checks their five original washout curves. These fixtures are not used as expected values for the final 5D/E data. The independent verifier recomputes all 18 current 4A trajectories directly from node coordinates, checks mean/sample SD and plot means, validates all bud counts/statistics, and checks export dimensions, colors and 4A/4B frame alignment. Historical 4A curves remain archived, not substituted for current final data.
