FIGURE 5C: WASHOUT BUD COUNTS

Run from the repository root after installing requirements.txt:

  python BB-94-washout/plotting_scripts/prepare_washout_data.py
  python BB-94-washout/plotting_scripts/plot_combined_log2_ratios_uniform.py
  python BB-94-washout/plotting_scripts/plot_combined_log2_ratios_symbols_no_legend.py

All defaults resolve relative to the scripts, so these commands also work from
another working directory when the script paths are supplied correctly. Repeated
runs replace generated results, while source files remain unchanged. Arial is the
original figure font. If unavailable, add --font "DejaVu Sans" to either plotting
command (font substitution can change text dimensions but not the data or tests).
Typography and line weights are shared through the root plot_style.py. Full
exported canvas sizes come from figure_sizes.json, with compact outside margins.
PDF and SVG are both produced by default, alongside a PNG preview.
Both marker variants are 1.00 x 1.40 inches including labels. Their BB-94 red is
firebrick (#B22222), matching Figure 5F; the color is shared in plot_palettes.json.
Control points use darkgray (#A9A9A9), matching Figure 5F, with black mean/SD
summary marks. Treatment/interval labels are tilted 45 degrees to provide a taller
plotting frame. X limits of -0.65 to 3.65 leave extra space between the outer
categories (at 0 and 3) and the frame. d1/d0 precedes
washout; d2/d1 is the subsequent washout interval. Fonts remain 7 pt for axis
labels and 6 pt for ticks, with smaller points to fit the gland-level observations.

The preparation command reads source_data/washout-inputs.json, checks the SHA-256
of every raw input, decodes image names, matches each gland across its three time
points, and calculates its two consecutive ratios and log2 ratios. A fingerprint
of the gland identities, times, and integer counts is compared with the validated
original figure data. It writes decoded and interval tables to results/<dataset>.

The manifest pairs these original, unmodified blinded-count inputs:

  210325-washout: 48 image counts, 16 glands (8 per treatment)
    210325-27-SMG-BB94-washout-1791559598-fileNameRecord.txt
    210325-27-SMG-BB94-washout-1791559598-renamed-8276-budCount.txt

  210401-washout: 78 image counts, 26 glands (13 per treatment)
    210401-03-SMG-BB94-washout-1791565130-fileNameRecord.txt
    210401-03-SMG-BB94-washout-1791565130-renamed-7422-budCount.txt

Each file is in its named dataset subfolder under source_data. The filename-record
files map original_file_name to converted_file_name; the count files supply
file_name and bud_number. These are the source measurements needed to reproduce
the plots; original microscopy movies are not needed for this analysis.
source-manifest-fig5.json records file sizes, SHA-256 hashes, and original
project-relative locations.

The first dataset is sampled at 0, 22, and 46 hours. The second uses day0/day1/day2;
these are stored as nominal 0/24/48 hours for pairing, without inferring precise
elapsed time. Its plate identifiers are part of gland identity. Datasets are
aligned by stage (d1/d0 and d2/d1) with no rate normalization. The pooled plot has
42 unique glands, 21 per treatment, and 84 interval values. Different treatments
are independent groups. A gland contributes only once to each interval test.

Both plot variants use two-sided Welch t-tests on per-gland log2 ratios with
Holm correction across the two control-versus-BB-94 interval comparisons.
Welch-Satterthwaite degrees of freedom allow unequal group variances. Tests
compare mean log2 ratios, assuming independent treatment groups and approximately
normal errors or an adequate mean approximation. Dataset and plate are not
covariates. The same gland contributes once to each interval test; the intervals
are not treated as four independent groups. Holm correction remains valid for
this correlated family. Mean and SD summarize the observations. The annotations
remain **** for both intervals. Adjusted P values are 8.419643671e-10 and
8.986676007e-10. Reported 95% mean-difference intervals are pointwise, not
multiplicity-adjusted. Error bar is SD.

Outputs in results/combined start with the current plot ID:

  Fig5C-washout-log2-ratios-uniform-stars.pdf/.svg/.png
  Fig5C-washout-log2-ratios-symbols-stars-no-legend.pdf/.svg/.png

Uniform uses 1.5 pt hollow circles with white interiors and 0.35 pt colored
outlines, matching the circle style in 4E/F/G; mean/SD marks are drawn above
the points. The other version retains filled circles for 210325-27 and filled
triangles for 210401-03, without a dataset legend. Both use gray for
control and firebrick red for BB-94. Control mean/SD marks are black; BB-94
mean/SD marks remain gray. The symbol explanation belongs in the figure legend.
All PDF and SVG text is converted to vector outlines, with no fonts needed to
view or place exports. Edit label wording in the scripts and regenerate the
files; exported glyphs are editable as shapes. Each variant includes a readable -statistics.txt
report, -statistics.csv test table, and -record.json method/provenance record.
Shared -data.csv, -by-dataset.csv, and -summary.csv tables retain every plotted
point and both dataset-specific and pooled summaries.

Each plotting command accepts --output-dir to relocate its outputs. The
preparation command accepts --output-root, but plotting defaults always read
results/<dataset>; standard reproduction should use the defaults above.
prepare_washout_data.py does not run the older per-dataset t tests. The retained
bud_count_stats.py module only matches glands and calculates interval ratios.
The older notebook and obsolete plot/statistics variants are intentionally not
part of this reproduction package.

Checks (from the repository root):

  python -m unittest discover -s BB-94-washout/plotting_scripts/tests -v

These verify shuffled identity-based decoding, timing aliases, plate identity,
invalid/duplicate inputs, missing observations, Welch standard errors and degrees
of freedom, Holm correction, and agreement with the validated figure inputs.
