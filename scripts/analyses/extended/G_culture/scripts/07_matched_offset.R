# Sensitivity (new, not in the manuscript): isolate-minus-MAG offset with phylogeny held constant.
# (a) within GTDB species containing both isolates and MAGs; (b) each MAG against its nearest isolate
# on the bac120 tree (patristic). Corrected labels; RGM and SGM.
suppressMessages(library(ape))
G <- "/path/to/project/analyses/extended/G_culture/"
P <- "/path/to/project/scripts/inputs/"
s <- read.delim(paste0(P, "tsv/S1_genome_metadata.tsv"))
cc <- read.delim(paste0(G, "corrected_culture_status.tsv")); s$new <- cc$new_status[match(s$tip, cc$tip)]
s$c3 <- sub(" .*", "", s$clade); s$iso <- grepl("^isolate", s$new)
s$iso_old <- s$culture_status %in% c("Cultured isolate", "NTM (isolate)", "Pathogen (isolate)")
tr <- read.tree(paste0(P, "myco_bac120.rooted.treefile"))
D <- cophenetic(tr)
out <- list()
for (lab in c("old", "new")) for (cl in c("RGM", "SGM")) {
  x <- s[s$c3 == cl, ]; isov <- if (lab == "new") x$iso else x$iso_old
  for (col in c("genome_size_Mb", "GC_pct", "coding_density_pct_uniform_prodigal")) {
    sp <- split(seq_len(nrow(x)), x$species)
    dsp <- unlist(lapply(sp, function(i) { a <- x[[col]][i][isov[i]]; b <- x[[col]][i][!isov[i]]
      if (sum(!is.na(a)) && sum(!is.na(b))) median(a, na.rm = TRUE) - median(b, na.rm = TRUE) else NULL }))
    I <- x$tip[isov]; M <- x$tip[!isov]
    nn <- I[apply(D[M, I, drop = FALSE], 1, which.min)]
    dn <- x[[col]][match(nn, x$tip)] - x[[col]][match(M, x$tip)]
    spire <- x$origin[match(M, x$tip)] == "SPIRE MAG"
    out[[length(out) + 1]] <- data.frame(labels = lab, clade = cl, metric = col,
      species_with_both = length(dsp), within_species_median_offset = round(median(dsp), 2),
      nn_pairs_allMAG = sum(!is.na(dn)), nn_median_offset_allMAG = round(median(dn, na.rm = TRUE), 2),
      nn_pairs_SPIRE = sum(!is.na(dn[spire])), nn_median_offset_SPIRE = round(median(dn[spire], na.rm = TRUE), 2),
      nn_median_patristic = round(median(apply(D[M, I, drop = FALSE], 1, min)), 3))
  }
}
o <- do.call(rbind, out); print(o, row.names = FALSE)
write.table(o, paste0(G, "S3_matched_offset_sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
