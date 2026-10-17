#!/usr/bin/env Rscript
# Table S17: phylogenetically corrected tests of FSM membership on the FULL
# 824-tip ingroup (FSM against RGM + SGM). Continuous traits use phylolm under a
# Pagel's lambda model; the binary opsin trait uses phyloglm (logistic_MPLE).
# Coding density is the uniform Prodigal call, the column used everywhere in the
# manuscript. Companion model on the restricted 525-tip FSM+RGM tree: Table S23.
suppressPackageStartupMessages({library(ape); library(phylolm)})
set.seed(1)
tr <- read.tree("inputs/myco_bac120.rooted.treefile")
d  <- read.delim("inputs/tsv/S1_genome_metadata.tsv",
                 stringsAsFactors = FALSE, na.strings = c("NA","","NaN","None"))
d$cl <- ifelse(grepl("^FSM", d$clade), "FSM",
        ifelse(grepl("^RGM", d$clade), "RGM",
        ifelse(d$clade == "SGM", "SGM", NA)))
d <- d[!is.na(d$cl) & d$tip %in% tr$tip.label, ]
tt <- drop.tip(tr, setdiff(tr$tip.label, d$tip))
rownames(d) <- d$tip; d <- d[tt$tip.label, ]
d$FSM <- as.integer(d$cl == "FSM")
cat(sprintf("ingroup tips = %d (FSM %d, RGM %d, SGM %d)\n",
            nrow(d), sum(d$cl=="FSM"), sum(d$cl=="RGM"), sum(d$cl=="SGM")))

cont <- list("genome size (Mb)"        = "genome_size_Mb",
             "coding density (%)"      = "coding_density_pct_uniform_prodigal",
             "acidic residues (%)"     = "acidic_residues_pct",
             "GC (%)"                  = "GC_pct")
res <- data.frame()
for (nm in names(cont)) {
  v <- cont[[nm]]
  dd <- d[!is.na(d[[v]]), ]
  tt2 <- drop.tip(tt, setdiff(tt$tip.label, dd$tip))
  dd <- dd[tt2$tip.label, ]
  y <- as.numeric(dd[[v]])
  df <- data.frame(y = y, FSM = dd$FSM); rownames(df) <- dd$tip
  fit <- phylolm(y ~ FSM, data = df, phy = tt2, model = "lambda")
  s <- summary(fit)$coefficients
  cat(sprintf("%-22s n=%3d  coef %+8.3f  lambda %.3f  p %.4g\n",
              nm, nrow(dd), s["FSM","Estimate"], fit$optpar, s["FSM","p.value"]))
  res <- rbind(res, data.frame(trait=nm, n=nrow(dd), coef=s["FSM","Estimate"],
                               lambda=fit$optpar, p=s["FSM","p.value"]))
}
dd <- d[!is.na(d$rhodopsin), ]
tt2 <- drop.tip(tt, setdiff(tt$tip.label, dd$tip)); dd <- dd[tt2$tip.label, ]
yb <- as.integer(dd$rhodopsin == "present")
gdf <- data.frame(y = yb, FSM = dd$FSM); rownames(gdf) <- dd$tip
gfit <- phyloglm(y ~ FSM, data = gdf, phy = tt2, method = "logistic_MPLE", btol = 30)
gs <- summary(gfit)$coefficients
cat(sprintf("%-22s n=%3d  beta %+8.3f  alpha %.3f  p %.4g\n",
            "opsin presence", nrow(dd), gs["FSM","Estimate"], gfit$alpha, gs["FSM","p.value"]))
res <- rbind(res, data.frame(trait="opsin presence", n=nrow(dd),
                             coef=gs["FSM","Estimate"], lambda=NA, p=gs["FSM","p.value"]))
dir.create("run_output", showWarnings=FALSE)
write.table(res, "run_output/S8_ingroup_models.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)
