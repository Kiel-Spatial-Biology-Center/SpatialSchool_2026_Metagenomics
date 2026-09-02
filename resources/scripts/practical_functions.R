# ==============================================================================
# Helper functions for the micro-scale spatial metagenomics practical
#
# Every function here is called from one of the .Rmd chapters. They are kept
# short and single-purpose so you can open this file and read what a step
# actually does.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. From mapping output to a community table
# ------------------------------------------------------------------------------

# breadth_of_coverage()
# Fraction of each genome's length that received at least one read, per
# microsample. A genome can collect reads simply because its sequence resembles
# something that is genuinely there, and those spurious reads land scattered
# across the genome. Breadth separates the two cases: a genome that is really
# present is covered evenly, a genome that only catches stray reads is not.
breadth_of_coverage <- function(covered_bases, genome_metadata) {
  stopifnot(identical(covered_bases$genome, genome_metadata$genome))
  covered_bases %>%
    mutate(across(where(is.numeric), ~ .x / genome_metadata$length))
}


# filter_by_coverage()
# Set a genome's read count to zero in any microsample where less than
# `min_coverage` of its length is covered. This is the single most important
# quality step for micro-scale data: microsamples hold very little DNA, so
# without it the rare tail of the community is mostly mapping noise.
filter_by_coverage <- function(read_counts, coverage, min_coverage = 0.3) {
  stopifnot(identical(read_counts$genome, coverage$genome))
  coverage %>%
    mutate(across(where(is.numeric), ~ ifelse(.x > min_coverage, 1, 0))) %>%
    mutate(across(-1, ~ .x * read_counts[[cur_column()]]))
}


# to_genome_counts()
# Convert read counts into genome counts (how many genome copies the reads
# imply). Long genomes capture more reads than short ones at equal abundance,
# so dividing by genome length in read units removes that bias.
to_genome_counts <- function(read_counts, genome_metadata, read_length = 150) {
  stopifnot(identical(read_counts$genome, genome_metadata$genome))
  read_counts %>%
    mutate(across(where(is.numeric), ~ .x / (genome_metadata$length / read_length)))
}


# alpha_diversity()
# Hill numbers per microsample: q = 0 is richness (counts genomes), q = 1 is
# neutral diversity (weights them by abundance), and the phylogenetic version
# of q = 1 additionally discounts genomes that are close relatives.
# Caution: with unnormalised pooling these correlate strongly with sequencing
# depth (rho > 0.9 for richness), partly because deeper sequencing detects more
# and partly because more biomass really does hold more genomes. The two cannot
# be separated in this data.
alpha_diversity <- function(genome_counts, tree) {
  counts <- genome_counts %>% column_to_rownames("genome")
  metrics <- list(
    richness     = hilldiv(counts, q = 0),
    neutral      = hilldiv(counts, q = 1),
    phylogenetic = hilldiv(counts, q = 1, tree = tree)
  )
  lapply(names(metrics), function(m) {
    metrics[[m]] %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("microsample") %>%
      rename(!!sym(m) := 2)
  }) %>%
    reduce(full_join, by = "microsample")
}


# ------------------------------------------------------------------------------
# 2. Preparing one cryosection for spatial analysis
# ------------------------------------------------------------------------------

# prepare_section()
# Pull out one cryosection and replace its zeros. Genome counts are
# compositional (only ratios are meaningful) and the log-ratio transforms used
# later are undefined at zero, so zeros are replaced by small positive values
# with a Bayesian-multiplicative estimator. Genomes too sparse to estimate
# (more than `z_warning` zeros) are dropped, which is why the returned table is
# narrower than the one you put in.
prepare_section <- function(section_id, genome_counts, metadata, z_warning = 0.95) {

  section_samples <- metadata %>%
    filter(cryosection == section_id) %>%
    pull(microsample)

  comm <- genome_counts %>%
    column_to_rownames("genome") %>%
    select(any_of(section_samples)) %>%
    t() %>%
    as.data.frame()

  # microsamples and genomes that are entirely empty carry no information
  comm <- comm[rowSums(comm) > 0, colSums(comm) > 0, drop = FALSE]

  comm_zero_replaced <- cmultRepl(comm, method = "GBM", output = "prop",
                                  z.warning = z_warning)

  comm_red <- comm[rownames(comm_zero_replaced), colnames(comm_zero_replaced),
                   drop = FALSE]

  metadata_section <- metadata %>%
    filter(microsample %in% rownames(comm_red)) %>%
    arrange(match(microsample, rownames(comm_red)))

  stopifnot(identical(metadata_section$microsample, rownames(comm_red)))

  list(
    section        = section_id,
    comm           = comm,       # raw genome counts, microsamples x genomes
    comm_red       = comm_red,   # same, restricted to genomes that survived
    comm_zero_repl = comm_zero_replaced,
    metadata       = metadata_section
  )
}


# neighbour_graph()
# Turn microsample coordinates into a neighbour network. A Gabriel graph joins
# two microsamples when no third one sits inside the circle that has them as
# its diameter, which gives a sensible "who is next to whom" for irregularly
# placed points without having to pick a distance cutoff.
neighbour_graph <- function(coords) {
  graph2nb(gabrielneigh(as.matrix(coords)), sym = TRUE)
}


# moran_correlogram()
# Moran's I measures whether values at neighbouring microsamples resemble each
# other more than expected by chance. Computing it at increasing lag orders
# (direct neighbours, neighbours-of-neighbours, ...) shows how far the
# resemblance reaches. Positive I at short lags = spatial structure.
moran_correlogram <- function(nb, values, order = 8) {
  cg <- sp.correlogram(nb, values, order = order, method = "I", zero.policy = TRUE)
  as.data.frame(cg$res) %>%
    setNames(c("morans_i", "expectation", "variance")) %>%
    mutate(
      lag   = row_number(),
      sd    = sqrt(variance),
      lower = morans_i - 1.96 * sd,
      upper = morans_i + 1.96 * sd
    )
}


# ------------------------------------------------------------------------------
# 3. Extended RLQ
# ------------------------------------------------------------------------------

# spatial_rlq()
# The core analysis. RLQ finds the axes along which the environment (R) and the
# species traits (Q) covary through the community table (L). The extended
# version (ESLTP) adds two more tables: the spatial arrangement of the
# microsamples and the phylogeny of the genomes. One axis therefore summarises,
# at once, where in the tissue a community sits, how diverse and deeply
# sequenced it is, and which kind of bacteria dominate it.
#
#   comm_red         genome counts for one cryosection (microsamples x genomes)
#   metadata_section its metadata, in the same row order
#   trait_pcoa       PCoA of genome functional traits (shared across sections)
#   tree             phylogeny of the genomes
#   root             exponent on the relative abundances; 0.5 = Hellinger
spatial_rlq <- function(section_id, comm_red, metadata_section, trait_pcoa, tree,
                        root = 0.5, correlogram_order = 8) {

  # --- L: the community table -------------------------------------------------
  # Relative abundances, then a square-root (Hellinger) transform that stops a
  # handful of dominant genomes from driving the whole ordination.
  comp <- decostand(comm_red, MARGIN = 1, method = "total")^root
  colnames(comp) <- gsub("\\.", ":", colnames(comp))

  # --- E: what we know about each microsample ---------------------------------
  # log_seq_depth is not a nuisance variable here. MSSM libraries are pooled
  # without equimolar normalisation, so the reads a microsample receives track
  # the DNA that went into it, and therefore its microbial biomass. Including it
  # lets RLQ give biomass an axis of its own instead of leaving it to
  # contaminate the axes we want to read spatially.
  env <- data.frame(
    log_seq_depth = log(metadata_section$after_filtering.total_bases),
    richness      = metadata_section$richness,
    log_host_dist = log(metadata_section$distance_host)
  )

  # --- T and P: what we know about each genome --------------------------------
  trait_scores <- trait_pcoa$vectors[
    rownames(trait_pcoa$vectors) %in% colnames(comp), 1:2, drop = FALSE]

  phy <- drop.tip(tree, setdiff(tree$tip.label, rownames(trait_scores)))

  # every table must list the genomes in the same order
  comp <- comp[, match(phy$tip.label, colnames(comp)), drop = FALSE]
  trait_scores <- as.data.frame(
    trait_scores[match(phy$tip.label, rownames(trait_scores)), , drop = FALSE])
  stopifnot(identical(phy$tip.label, colnames(comp)),
            identical(phy$tip.label, rownames(trait_scores)))

  phylog <- newick2phylog(write.tree(phy))
  colnames(comp) <- gsub(":", "_", colnames(comp))
  rownames(trait_scores) <- gsub(":", "_", rownames(trait_scores))

  # --- S: the spatial arrangement ---------------------------------------------
  spa <- metadata_section[, c("Xcoord", "Ycoord")]
  nb  <- neighbour_graph(spa)

  correlograms <- list(
    richness      = moran_correlogram(nb, log(env$richness), correlogram_order),
    log_seq_depth = moran_correlogram(nb, env$log_seq_depth, correlogram_order),
    log_host_dist = moran_correlogram(nb, env$log_host_dist, correlogram_order)
  )

  # --- the five ordinations RLQ combines --------------------------------------
  # dudiL first: it supplies the row and column weights that make the other
  # four ordinations comparable to it.
  coacomp <- dudi.coa(comp, scan = FALSE, nf = ncol(comp))

  spatial_vectors <- scores.neig(nb2neig(nb))
  pcaspa <- dudi.pca(spatial_vectors, row.w = coacomp$lw, scan = FALSE,
                     nf = ncol(spatial_vectors))
  pcaenv <- dudi.pca(env, row.w = coacomp$lw, scannf = FALSE, nf = 2)

  pcophy <- dudi.pco(
    as.dist(as.matrix(phylog$Wdist)[names(comp), names(comp)]),
    coacomp$cw, full = TRUE)

  trait_dist <- dist.ktab(ktab.list.df(list(trait_scores)), c("Q"), scan = FALSE)
  pcotraits  <- dudi.pco(trait_dist, coacomp$cw, full = TRUE)

  rlqmix <- rlqESLTP(pcaenv, pcaspa, coacomp, pcotraits, pcophy,
                     scan = FALSE, nf = 2)

  list(
    section      = section_id,
    comp         = comp,
    env          = env,
    tree         = phy,
    coords       = spa,
    nb           = nb,
    correlograms = correlograms,
    rlqmix       = rlqmix,
    eig_prop     = rlqmix$eig / sum(rlqmix$eig)
  )
}


# flip_rlq_axes()
# The sign of an ordination axis is arbitrary: flipping it changes nothing about
# the result, only which end of the gradient is drawn in red. We flip chosen
# axes so that the caecum and colon figures point the same way and can be
# compared side by side.
flip_rlq_axes <- function(rlq_result, axes) {
  slots <- c("lR", "mR", "lQ", "mQ",
             "lR_givenE", "lR_givenS", "lQ_givenT", "lQ_givenP")
  for (ax in axes) {
    for (s in slots) {
      m <- rlq_result$rlqmix[[s]]
      if (!is.null(m) && ncol(m) >= ax) rlq_result$rlqmix[[s]][, ax] <- -m[, ax]
    }
  }
  rlq_result
}


# rlq_site_scores()
# Where each microsample falls on the RLQ axes, joined to its metadata. RLQ
# labels its site rows by position, so we map them back onto microsample names
# through the community table. One row per microsample.
rlq_site_scores <- function(rlq_result, metadata_section) {
  lR <- as.data.frame(rlq_result$rlqmix$lR)
  microsamples <- rownames(rlq_result$comp)[as.numeric(rownames(lR))]

  tibble(microsample = microsamples) %>%
    bind_cols(lR %>% setNames(paste0("axis", seq_len(ncol(lR))))) %>%
    left_join(metadata_section, by = "microsample")
}


# rlq_genome_scores()
# Where each genome falls on the same axes, joined to its taxonomy and size.
# One row per genome. RLQ rewrites ":" as "_" in genome names, so we match them
# back against the catalogue rather than guessing at the substitution.
rlq_genome_scores <- function(rlq_result, genome_metadata) {
  lQ <- as.data.frame(rlq_result$rlqmix$lQ)
  pool <- genome_metadata$genome
  genomes <- pool[match(rownames(lQ), gsub(":", "_", pool))]

  tibble(genome = genomes) %>%
    bind_cols(lQ %>% setNames(paste0("axis", seq_len(ncol(lQ))))) %>%
    left_join(genome_metadata, by = "genome")
}


# ------------------------------------------------------------------------------
# 4. Putting results back on the tissue
# ------------------------------------------------------------------------------

# cryosection_image()
# Load the brightfield photograph of one cryosection.
cryosection_image <- function(section_id, dir = "resources/images", brightness = 100) {
  path <- file.path(dir, paste0(section_id, "_bright.png"))
  if (!file.exists(path)) return(NULL)
  magick::image_modulate(magick::image_read(path), brightness = brightness)
}


# plot_on_section()
# Draw the microsamples at their real positions on the cryosection photograph,
# coloured by any per-microsample value. This is the plot that makes a spatial
# gradient visible: an RLQ axis painted here should show structure, a shuffled
# version of the same values should not.
#
#   diverging = TRUE  red-white-blue scale centred on zero (RLQ axes)
#   diverging = FALSE viridis scale binned by quantile (richness, depth, ...)
plot_on_section <- function(section_id, values_df, value_col, value_label,
                            pixel_coords, diverging = TRUE, n_bins = 5,
                            image_size = 1000, dir = "resources/images") {

  img <- cryosection_image(section_id, dir = dir)
  if (is.null(img)) {
    message("no image for ", section_id)
    return(invisible(NULL))
  }

  # image pixel y counts down from the top, ggplot's y counts up from the bottom
  plot_df <- pixel_coords %>%
    filter(cryosection == section_id, !is.na(pixel_x), !is.na(pixel_y)) %>%
    inner_join(values_df %>% select(microsample, all_of(value_col)),
               by = "microsample") %>%
    mutate(pixel_y_flip = image_size - pixel_y)

  if (diverging) {
    lim <- max(abs(plot_df[[value_col]]), na.rm = TRUE)
    colour_aes <- aes(x = pixel_x, y = pixel_y_flip, colour = .data[[value_col]])
    colour_scale <- scale_colour_fermenter(
      palette = "RdBu", n.breaks = 7, limits = c(-lim, lim), name = value_label)
  } else {
    breaks <- unique(quantile(plot_df[[value_col]],
                              probs = seq(0, 1, length.out = n_bins + 1),
                              na.rm = TRUE))
    if (length(breaks) < 3) {
      breaks <- pretty(range(plot_df[[value_col]], na.rm = TRUE), n = n_bins)
    }
    plot_df$.bin <- cut(plot_df[[value_col]], breaks = breaks,
                        include.lowest = TRUE, dig.lab = 4)
    colour_aes <- aes(x = pixel_x, y = pixel_y_flip, colour = .bin)
    colour_scale <- scale_colour_viridis_d(name = value_label, option = "viridis")
  }

  ggplot(plot_df, colour_aes) +
    annotation_raster(as.raster(img), xmin = 0, xmax = image_size,
                      ymin = 0, ymax = image_size) +
    geom_point(size = 2) +
    colour_scale +
    coord_fixed(xlim = c(0, image_size), ylim = c(0, image_size), expand = FALSE) +
    labs(title = section_id) +
    theme_void() +
    theme(plot.title = element_text(size = 10, hjust = 0.5))
}


# plot_neighbour_graph()
# The Gabriel graph drawn over the microsample coordinates: the network that
# Moran's I and the spatial component of RLQ actually use.
plot_neighbour_graph <- function(coords, nb, point_values = NULL,
                                 value_label = NULL) {
  coords <- as.data.frame(coords) %>% setNames(c("x", "y"))
  edges <- do.call(rbind, lapply(seq_along(nb), function(i) {
    j <- nb[[i]]
    j <- j[j > i]                       # each edge once
    if (length(j) == 0 || identical(j, 0L)) return(NULL)
    data.frame(x = coords$x[i], y = coords$y[i],
               xend = coords$x[j], yend = coords$y[j])
  }))

  p <- ggplot() +
    geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend),
                 colour = "grey70", linewidth = 0.3)

  if (is.null(point_values)) {
    p <- p + geom_point(data = coords, aes(x, y), size = 1.6, colour = "grey20")
  } else {
    p <- p +
      geom_point(data = coords %>% mutate(value = point_values),
                 aes(x, y, colour = value), size = 1.8) +
      scale_colour_viridis_c(name = value_label)
  }

  p + coord_fixed() + theme_minimal() +
    labs(x = "X (µm)", y = "Y (µm)")
}


# plot_tree_gradient()
# The genome phylogeny with tips coloured by their RLQ score. If the axis
# separates whole clades, the gradient is phylogenetically conserved; if the
# colours are scattered across the tree, it is not.
plot_tree_gradient <- function(tree, genome_scores, axis_col, axis_label,
                               tip_colour_by = "order", colours = NULL) {

  tree_data <- genome_scores %>%
    mutate(label = genome) %>%
    select(label, all_of(axis_col), all_of(tip_colour_by), phylum, family, length)

  lim <- max(abs(tree_data[[axis_col]]), na.rm = TRUE)

  p <- ggtree(tree, ladderize = FALSE) %<+% tree_data +
    geom_tippoint(aes(fill = .data[[axis_col]]), shape = 21, size = 2.4,
                  stroke = 0.2, colour = "grey30") +
    scale_fill_fermenter(palette = "RdBu", n.breaks = 7,
                         limits = c(-lim, lim), name = axis_label) +
    theme(legend.position = "right")

  if (!is.null(colours)) {
    p <- p + geom_tiplab(aes(label = .data[[tip_colour_by]],
                             colour = .data[[tip_colour_by]]),
                         size = 1.8, offset = 0.02) +
      scale_colour_manual(values = colours, guide = "none")
  }
  p
}
