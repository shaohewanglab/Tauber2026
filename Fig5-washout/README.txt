FIGURE 5: WASHOUT BUD COUNTS

Run from the repository root after installing requirements.txt:

  python Fig5-washout/plotting_scripts/prepare_washout_data.py
  python Fig5-washout/plotting_scripts/plot_combined_log2_ratios_uniform.py
  python Fig5-washout/plotting_scripts/plot_combined_log2_ratios_symbols_no_legend.py

All defaults resolve relative to the scripts, so these commands also work from
another working directory when the script paths are supplied correctly. Repeated
runs replace generated results, while source files remain unchanged. Arial is the
original figure font. If unavailable, add --font "DejaVu Sans" to either plotting
command (font substitution can change text dimensions but not the data or tests).

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

Both plot variants use two-sided Mann-Whitney U tests with Holm correction across
the two control-versus-BB-94 interval comparisons. SciPy method="auto" and
use_continuity=True reproduce the agreed test; both current comparisons use the
asymptotic method with tied ranks and continuity correction. Dataset and plate are
not covariates. Mean and SD are descriptive summaries; Mann-Whitney tests compare
distributions. The annotations remain **** for both intervals. Error bar is SD.

Outputs in results/combined retain the current filenames:

  combined-washout-log2-ratios-uniform-stars.pdf/.svg/.png
  combined-washout-log2-ratios-symbols-stars-no-legend.pdf/.svg/.png

Uniform uses circles for every gland. The other version uses circles for
210325-27 and triangles for 210401-03, without a dataset legend. Both use black for
control and red for BB-94. The symbol explanation belongs in the figure legend.
SVG text remains editable. Each variant includes a readable -statistics.txt
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

  python -m unittest discover -s Fig5-washout/plotting_scripts/tests -v

These verify shuffled identity-based decoding, timing aliases, plate identity,
invalid/duplicate inputs, missing observations, exact and tied Mann-Whitney
calculations, Holm correction, and agreement with the validated figure inputs.
