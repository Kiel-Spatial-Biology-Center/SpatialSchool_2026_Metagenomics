# Micro-scale spatial metagenomics — student practical

Practical for the **Kiel Spatial Biology Summer Course 2026**.

A teaching version of the analysis behind *"Micro-scale spatial metagenomics:
revealing high-resolution spatial biogeography of gut microbiomes"*
([3D'omics](http://www.3domics.eu), Horizon 2020).

**Read the book online:**
<https://kiel-spatial-biology-center.github.io/SpatialSchool_2026_Metagenomics>

It keeps only what is needed to understand **how spatial structure in a gut
microbiome is detected and tested**, and drops the method validation,
benchmarking and strain-level analyses of the published work.

## Contents

| Chapter | Question | Methods |
|---|---|---|
| 1 · Data | What does a microsample community look like? | Breadth-of-coverage filtering, genome counts, Hill numbers |
| 2 · Space | Is the community spatially structured? | Tissue maps, CLR, PCA, PERMANOVA, Gabriel graphs, Moran's I |
| 3 · RLQ | What explains the structure? | Extended RLQ (environment + space + traits + phylogeny) |
| 4 · Exercises | Does it hold up? | Thresholds, null models, table ablation, interpretation |

Chapter 2 asks whether the community is spatially structured, using the standard
microbial ecology toolkit. Chapter 3 comes back to the same data with extended
RLQ, which maximises covariance between what we know about microsamples and what
we know about genomes. Chapter 4 tests how much of the answer survives changed
thresholds, null models and dropped tables.

The book also takes seriously a protocol detail that is unusual for
metagenomics: MSSM libraries are **not pooled equimolarly**. Equal microsample
areas get equal PCR cycles and are pooled as they come, with one cryosection per
sequencing run. Sequencing depth therefore tracks input biomass instead of being
equalised away, and composition closes at the cryosection level. This makes
depth a quasi-quantitative biological variable rather than a nuisance. Chapter 1
verifies the design in the data; Chapters 2 and 3 build on it.

## Data

Everything is in this repository — **nothing to download**.

| File | Size | Contents |
|---|---|---|
| `resources/data/mssm_practical_data.Rdata` | 0.4 MB | 336 microsamples × 223 genomes: counts, coverage, metadata, taxonomy, phylogeny, metabolic traits |
| `resources/images/` | 1.7 MB | Two cryosection photographs + microsample pixel coordinates |

Two cryosections from the same chicken (G121), so region differences are not
confounded with individual:

| Cryosection | Region | Microsamples |
|---|---|---|
| `G121eI117A` | Caecum | 168 |
| `G121eO306A` | Colon | 168 |

The 223 genomes are MAGs assembled beforehand from conventional macro-scale
sequencing of the same animals; microsample reads are mapped against that
catalogue.

## Running it

The rendered HTML book is available at
<https://kiel-spatial-biology-center.github.io/SpatialSchool_2026_Metagenomics>,
so it can be followed without running anything. To build it locally, open
`MSSM_practical.Rproj` in RStudio and knit, or:

```r
bookdown::render_book(input = ".", output_format = "bookdown::gitbook",
                      output_dir = "docs")
```

The whole book renders in a couple of minutes; each RLQ takes about a second, so
you can re-run and experiment freely.

Chapters 2 and 3 rebuild what they need from earlier chapters, so each can also
be run on its own after `index.Rmd`.

### Requirements

R ≥ 4.2. Packages install on first run via `pacman` (see `index.Rmd`):

- **CRAN** — tidyverse, ape, phytools, vegan, zCompositions, ade4, adespatial,
  adegraphics, spdep, patchwork, ggrepel, magick
- **Bioconductor** — ggtree
- **GitHub** — [hilldiv2](https://github.com/anttonalberdi/hilldiv2),
  [distillR](https://github.com/anttonalberdi/distillR)

## Layout

```
index.Rmd            introduction, setup, palettes, data loading
01-data.Rmd          from mapping output to a community table
02-space.Rmd         testing for spatial structure (and failing)
03-rlq.Rmd           extended RLQ (and succeeding)
04-exercises.Rmd     exercises, in six sets
resources/
  data/              the trimmed dataset
  images/            cryosection photographs and pixel coordinates
  scripts/
    practical_functions.R   all helpers, commented for reading
    rlqESLTP.R              extended RLQ (Pavoine et al. 2011)
```

`resources/scripts/practical_functions.R` is meant to be read. Each function
carries a comment explaining *why* the step exists, not just what it does.

## What was left out

Relative to the published analysis: microsample-area validation and detection
limits, lysis-condition and resource-optimisation experiments, read-level QC
pipelines (fastp, Kraken, SingleM, Nonpareil), host and environmental
contamination assessment, macro- vs micro-scale MAG catalogue comparison,
differential abundance, FISH and macro-scale metagenomics cross-validation, mock
communities, pangenome and SNP-level strain analyses, and the six-cryosection
replication.

Full analysis: [github.com/3d-omics/MSSM](https://github.com/3d-omics/MSSM)

## Citation

Pietroni C., Langa J., Odriozola I., Bogri A., Alberdi A. *Micro-scale spatial
metagenomics: revealing high-resolution spatial biogeography of gut
microbiomes.* https://doi.org/10.1101/2025.09.30.679663

The manuscript is currently under review in PNAS; the preprint above is the
version to cite in the meantime.
