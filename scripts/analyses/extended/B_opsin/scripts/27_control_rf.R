#!/usr/bin/env Rscript
# nRF to the species tree for opsin and control-gene ML trees on the same genome sets (common taxa only)
suppressPackageStartupMessages({library(ape); library(phangorn)})
O <- "/path/to/project/analyses/extended/B_opsin/reconciliation"
sp <- read.tree("/path/to/project/scripts/inputs/myco_bac120.rooted.treefile"); sp$node.label <- NULL
rows <- list()
for (s in c("A_all","A_FSM","B_all","B_FSM")) {
  op <- read.tree(file.path(O, sprintf("au_%s_ml.treefile", s)))
  th <- read.tree(file.path(O, sprintf("control/thrS_%s_ml.treefile", s))); tp <- read.tree(file.path(O, sprintf("control/topA_%s_ml.treefile", s)))
  com <- Reduce(intersect, list(op$tip.label, th$tip.label, tp$tip.label))
  f <- function(g) RF.dist(unroot(keep.tip(g, com)), unroot(keep.tip(sp, com)), normalize=TRUE)
  rows[[s]] <- data.frame(set=s, common_genomes=length(com), nRF_opsin=round(f(op),3), nRF_thrS=round(f(th),3), nRF_topA=round(f(tp),3))
}
R <- do.call(rbind, rows); print(R, row.names=FALSE)
write.table(R, file.path(O, "control/control_rf.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
