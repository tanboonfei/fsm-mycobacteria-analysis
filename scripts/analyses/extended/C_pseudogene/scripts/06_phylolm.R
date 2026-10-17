# Phylogenetic regression (Pagel's lambda, as in Table S4f) of decay-screen metrics
# on FSM membership, SPIRE MAGs with completeness >= 95%, bac120 tree.
suppressPackageStartupMessages({library(ape); library(phylolm)})
D <- "/path/to/project/analyses/extended/C_pseudogene"
tr <- read.tree("/path/to/project/scripts/inputs/myco_bac120.rooted.treefile")
g <- read.delim(file.path(D, "results/per_genome_decay_screen_grouped.tsv"), na.strings = c("NA", "", "NaN", "None"))
run <- function(sub, label) {
  sub$tip <- ifelse(grepl("^spire", sub$genome), paste0(sub$genome, ".fa"), sub$tip)
  sub <- sub[sub$tip %in% tr$tip.label, ]; rownames(sub) <- sub$tip
  tt <- keep.tip(tr, sub$tip)
  sub$FSM <- as.integer(sub$clade == "FSM"); sub$lcpm <- log10(sub$contigs_per_Mb)
  out <- NULL
  for (y in c("pct_fs_or_stop", "pct_frameshift", "pct_internal_stop", "pct_orf_lt80", "pct_split", "pct_any_decay", "spacer_median_bp", "pct_spacer_gt300")) {
    for (f in c("FSM", "FSM + lcpm + completeness")) {
      fit <- phylolm(as.formula(paste(y, "~", f)), data = sub, phy = tt, model = "lambda")
      cf <- summary(fit)$coefficients
      out <- rbind(out, data.frame(panel = label, n_FSM = sum(sub$FSM), n_RGM = sum(1 - sub$FSM), trait = y, model = paste(y, "~", f),
        beta_FSM = cf["FSM", "Estimate"], SE = cf["FSM", "StdErr"], p = cf["FSM", "p.value"], lambda = fit$optpar))
    }
  }
  out
}
mag <- g[g$panel == "SPIRE_MAG" & g$completeness >= 95 & g$clade %in% c("FSM", "RGM"), ]
mt <- g[g$group %in% c("FSM SPIRE MAG (>=95%)", "RGM SPIRE MAG (>=95%), matched"), ]
res <- rbind(run(mag, "SPIRE MAGs >=95% complete, FSM vs all RGM"), run(mt, "SPIRE MAGs >=95%, FSM vs matched RGM"))
write.table(res, file.path(D, "results/phylolm.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(res, digits = 3)
