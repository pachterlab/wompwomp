# paper_figures

Scripts that regenerate the figures and reported objective values for
*"Optimizing alluvial plots"* (Rich, Oakes, Pachter).

These are the canonical figure-generation sources. The manuscript repo
(`wompwomp_paper`) holds **copies** of the finished images in
`figures_and_tables/`; regenerate here, then copy over.

## Contents

| script | paper figures |
| --- | --- |
| `examples-with-small-data.Rmd` | Fig. 1A–C (tissue → cluster), Supp. Table 1 |
| `examples-with-large-data.Rmd` | Fig. 1D–E (Game of Thrones), Fig. 3D–F (13-method clustering), Supp. Fig. 3–4, Supp. Fig. 6 |
| `examples-with-labeling.Rmd` | Fig. 2, Supp. Fig. 2 |
| `8-cubed.Rmd` | Fig. 3A–B (scRNA, two genes) |
| `8_cubed_sankey.ipynb` | Fig. 3C (Sankey) |
| `kits-vs-inst.Rmd` | Supp. Fig. 5 (KiTS19 model comparison) |
| `global-optimum-comparison.Rmd` | empirical checks in the text |
| `biowomp_intro.Rmd` | tutorial (not a paper figure) |

Each `.Rmd` writes its output to `paper_figures/output/`.

## Requirements

Use the Docker image built from `../Dockerfile` — it pins the toolchain
(libpng/freetype/harfbuzz, `ggfittext`, `ggrastr`, Bioconductor
`DuoClustering2018`) that the plotting and clustering-dataset steps need.

```
cd ..                 # wompwomp repo root
docker build -t wompwomp:figures .
docker run --rm -it -p 8787:8787 -e PASSWORD=changeme wompwomp:figures
# open http://localhost:8787 (user: rstudio), then knit the scripts in paper_figures/
```

`plot_alluvial()` and everything else these scripts need is in `wompwomp`
itself (the former `biowomp` package was merged in).

## Data

Two datasets are not downloadable and must be placed in `paper_figures/output/`
before running the corresponding script:

- `Akr1c21_Slc7a12_grouped.csv`  — for `8-cubed.Rmd`
- `kits_alluvial_data.csv`       — for `kits-vs-inst.Rmd`

Game of Thrones data downloads automatically; the 13-method clustering data
comes from `DuoClustering2018`.

## After regenerating

Copy the finished figure files from `paper_figures/output/` into
`wompwomp_paper/figures_and_tables/`, and update any objective values quoted in
`main.tex` that changed.
