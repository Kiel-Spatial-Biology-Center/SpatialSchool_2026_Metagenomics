suppressMessages({library(tidyverse); library(ape); library(vegan); library(zCompositions); library(ade4); library(adespatial); library(adegraphics); library(spdep); library(patchwork); library(magick); library(hilldiv2)})
select <- dplyr::select
setwd("/Users/anttonalberdi/Downloads/MSSM-practical")
theme_set(theme_classic(base_size = 12))
source("resources/scripts/rlqESLTP.R"); source("resources/scripts/practical_functions.R")
load("resources/data/mssm_practical_data.Rdata")
pixel_coords <- read_csv("resources/images/microsample_pixel_coords.csv", show_col_types = FALSE)

gcov <- breadth_of_coverage(genome_covered_bases, genome_metadata)
genome_counts <- to_genome_counts(filter_by_coverage(genome_read_counts, gcov, 0.3), genome_metadata)
sample_metadata <- sample_metadata %>% left_join(alpha_diversity(genome_counts, genome_tree), by = "microsample")
caecum <- suppressWarnings(prepare_section("G121eI117A", genome_counts, sample_metadata))
colon  <- suppressWarnings(prepare_section("G121eO306A", genome_counts, sample_metadata))
gu <- union(colnames(caecum$comm_red), colnames(colon$comm_red))
gifts_used <- genome_gifts[rownames(genome_gifts) %in% gu, !grepl("^S0|^D04", colnames(genome_gifts)), drop = FALSE]
trait_pcoa <- ape::pcoa(dist(gifts_used))
rlq_caecum <- flip_rlq_axes(spatial_rlq("G121eI117A", caecum$comm_red, caecum$metadata, trait_pcoa, genome_tree), c(1, 2))
sites <- rlq_site_scores(rlq_caecum, caecum$metadata)

mk <- function(col, lab) {
  plot_on_section("G121eI117A", sites, col, lab, pixel_coords) +
    theme(legend.position = "bottom",
          legend.key.width = unit(0.9, "cm"), legend.key.height = unit(0.28, "cm"),
          legend.title = element_text(size = 11), legend.text = element_text(size = 8),
          plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
}
p <- (mk("axis2", "Axis 2  ·  spatial gradient") + labs(title = "Axis 2 — community gradient")) +
     (mk("axis1", "Axis 1  ·  biomass")          + labs(title = "Axis 1 — biomass"))
ggsave("/private/tmp/claude-501/-Users-anttonalberdi-Downloads-mssm-for-teaching/26a26b5e-15b0-4057-9a53-ad702a563b7b/scratchpad/slide_axes.png",
       p, width = 11, height = 6.2, dpi = 150, bg = "white")
cat("done\n")
