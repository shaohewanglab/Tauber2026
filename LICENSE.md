# Licensing

Copyright © 2026 Tauber-2026-modeling contributors, except where an existing
copyright notice identifies another holder.

| Material | License |
| --- | --- |
| Original code, notebooks, software configuration, and associated documentation; exception below | [MIT](LICENSE-CODE) |
| `3D-RCAN-training/apply.py`, adapted from 3D-RCAN | [CC BY-NC 4.0](3D-RCAN-training/LICENSE-APPLY), noncommercial use only |
| Source and test data, results, data manifests, and validation records | [CC BY 4.0](LICENSE-DATA) |

The modified `apply.py` retains the upstream copyright notice and records Shaohe
Wang's changes in its header. `3D-RCAN-training/voxel_size_utils.py` has its own
MIT notice, which is preserved.

The mesh-morphometry plotting script adapts the scripts shared by Justin.
Original script and supplied-PDF hashes are recorded in
`mesh-morphometry/validation/provenance.json`.

The model-morphometry script consolidates Justin's supplied acircularity
scripts, retained unchanged with hashes in
`model-morphometry/source-data-manifest.json`.

The eight scripts in `image-processing/` were shared by Justin Tauber,
Soft Math Lab (Harvard University), in collaboration with Shaohe Wang.
They retain their original author headers and are preserved unchanged;
`image-processing/provenance.json` records their supplied-file hashes and
the hash of the original README, retained as `README-Justin.md`.

See [LICENSE-DATA](LICENSE-DATA) for data scope, attribution, and links to the
full CC BY 4.0 legal code and summary. Third-party dependencies retain their
respective licenses.
