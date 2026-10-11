FIGURE 4: BB-94 TITRATION

From the package root:

  python Fig4-BB94-titration/plotting_scripts/plot_bud_count_titration.py --output-dir Fig4-BB94-titration/results/bud-count --overwrite
  python Fig4-BB94-titration/plotting_scripts/plot_collagen_iv_titration.py --output-dir Fig4-BB94-titration/results/collagen-IV --overwrite

The copied JSON manifests resolve all inputs relative to their own locations.
Source tables are in ../source_data. Both CLIs default to Mann-Whitney U with
Holm correction over the six comparisons against DMSO (0), separately per plot.
The plotting style and numerical methods match the accepted manuscript plots.

Outputs retain their current names: Fig4-bud-count-titration and
Fig4-collagen-IV-intensity-titration, each as PDF/SVG/PNG plus report/data files.
The -comparisons.csv and -versus-control.csv files each contain the six tested
contrasts. -glands.csv contains the individual gland endpoints; -summary.csv
contains n, mean and SD; -run.json records inputs, code hashes and settings.

Input formats:
  Bud filename records: original_file_name, converted_file_name.
  Bud counts: file_name, bud_number. Tab or comma delimiters are supported.
  Collagen CSV: file_name (gland), id (scan), channel and AUC.

Bud identity is cohort + group + gland number, paired across days before forming
ratios. The manifests map dose and day explicitly. Unknown/missing/duplicate
matches, incomplete pairs and invalid counts cause errors. Only filename-marked
placeholders and non-endpoint days are excluded. Actual single-bud glands stay.

Collagen scans are averaged within each gland before testing. The control
normalization denominator is 262.485293693283 AUC units (equal gland weights),
matching the accepted CLI analysis. It differs by about 0.1% from the old
notebook's scan-weighted divisor; all current points and tests are reproduced.

Useful options:
  --config PATH       Input manifest for another dataset.
  --output-dir DIR    Required output directory.
  --prefix NAME       Output filename stem.
  --formats pdf svg png
  --font Arial        Use --font "DejaVu Sans" if Arial is unavailable.
  --width 1.9 --height 1.75
  --overwrite         Replace this prefix's generated outputs.

The existing --test tukey and --test games-howell modes are retained in the
extracted CLI for compatibility. They are not used in reproduce_plots.py or the
included accepted results; unlike the current Mann-Whitney mode, their families
contain all pairwise dose comparisons. Their additional helper module is included.

The original Jupyter notebooks and their sw_ dependencies are not required.
