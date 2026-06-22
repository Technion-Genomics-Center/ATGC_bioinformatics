# `rmd/` — example-figure sources

Each subfolder here is a self-contained R Markdown that **generates the example
figures** shown on the bioinformatics services page (`../index.html`).

Rationale: figures rendered by real R packages (ggplot2, ComplexHeatmap,
Seurat, DESeq2, etc.) look more convincing and realistic than charts hand-built
in JavaScript/SVG. We keep the `.Rmd` sources here so every figure is
reproducible and easy to re-style later.

## Convention

```
rmd/
  <section>/            e.g. scrna/, rnaseq/, methylation/, 16s/
    <section>.Rmd       script that produces the figures (seeded, synthetic data)
    figs/               output figures (SVG preferred; PNG fallback)
```

- Use a fixed `set.seed()` so figures are reproducible.
- Synthetic/representative data only — no patient/project data in this repo.
- Prefer **SVG** output (`svglite`/`ggsave(".svg")`) so figures stay crisp and
  can be inlined into `index.html`.

## Render

```bash
"/c/Program Files/R/R-4.5.1/bin/Rscript.exe" -e 'rmarkdown::render("rmd/<section>/<section>.Rmd")'
```

(Pandoc is bundled with RStudio — if rendering from the CLI, set
`RSTUDIO_PANDOC` in the outer call; see the project notes.)
