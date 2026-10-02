# Source and validation boundaries

The 18 Si/Al calculations and three targeted unit tests were run using the
unmodified CP2K mainline revision
[`dd94896b48d429648fb3640799c9d198ff8de1e3`](https://github.com/cp2k/cp2k/tree/dd94896b48d429648fb3640799c9d198ff8de1e3),
frozen on 16 September 2026. `validation/build_metadata.json` records the
executable SHA-256, compiler flags, thread count, and embedded source revision.
Per-case metadata identify the executable used for every retained output.

The `validation_additional/` cases (12 Materials-Project structures across four
Monkhorst-Pack/MacDonald meshes plus a five-mode same-mesh comparison on the
shifted and Gamma-centred $4^3$ grids) were run with a separate local build at
CP2K revision `305e089a6f` (with CUDA and libxc), so their absolute reference
energies are recorded but not hard-checked against the published Si/Al values.
`validation_additional/structures.json` stores both the as-given Materials-Project
metric (`lattice_asgiven`) and the metric snapped to the spglib-detected point
group; the latter is what the generated inputs and all reported counts use, so
that symmetry-reduced and full-grid runs integrate the same cell.
`validation_additional/results_additional.json` and
`symmetry_modes_results{,_gamma}.json` hold the parsed per-case results. The
`additional_tables.*` and `symmetry_modes_tables*.tex` files are the table
sources. The 168 raw CP2K `.out` logs are retained compressed as
`validation_additional/outs.tar.xz` (xz, ~1.2% of the uncompressed size);
extract with `tar -xf outs.tar.xz` to restore the original per-case paths.

`validation/source_history.json` is a title-filtered inventory of relevant
and adjacent mainline changes through that freeze. It is a source-history
record, not evidence that every listed feature was rerun for this article.
The SI's source-module and capability tables provide the corresponding
inspection map. Public CP2K source at the pinned revision is the authoritative
implementation; source-file copies are intentionally not bundled here.

The later [Loewdin DFT+U k-point implementation](https://github.com/cp2k/cp2k/commit/617be45fc17049d0e3d662055d936478c44379d8)
entered mainline on 25 September 2026. Its regression inputs cover restricted
and unrestricted spin, standard diagonalization and direct OT, including
numerical force and stress checks. It was not used for the Si/Al calculations
in this archive. Its k-point path requires a full mesh or time-reversal-only
reduction and does not support atomic point-group reduction.

The in-process Wannier90 library interface discussed in the manuscript was
inspected as a separate development extension. No Wannier90-library numerical
result is claimed from this dataset. The Si/Al tests establish fixed-mesh
equivalence for the stated cells and one-thread layout, not general parallel
scaling, nonzero-force accuracy on distorted cells, or downstream-method
validation. The retained failed two-thread logs document one unresolved
diagnostic at the validation snapshot; they are not successful test cases.
