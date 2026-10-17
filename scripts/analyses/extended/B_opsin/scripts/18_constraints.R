#!/usr/bin/env Rscript
# Partial constraint trees for AU tests: (i) species-tree splits with UFBoot >= 95 only; (ii) FSM opsins of the lineage monophyletic
suppressPackageStartupMessages(library(ape))
O <- "/path/to/project/analyses/extended/B_opsin"
P <- "/path/to/project/scripts"
sp <- read.tree(file.path(P,"inputs/myco_bac120.rooted.treefile"))
s1 <- read.delim(file.path(P,"inputs/tsv/S1_genome_metadata.tsv"), stringsAsFactors=FALSE)
fsm <- s1$tip[grepl("^FSM", s1$clade)]
for (nm in c("A_all","A_FSM","B_all","B_FSM")) {
  tips <- names(read.FASTA(file.path(O, sprintf("reconciliation/%s.aln", nm)), type="AA"))
  s <- keep.tip(sp, tips)
  sup <- suppressWarnings(as.numeric(s$node.label)); sup[is.na(sup)] <- 100
  s$edge.length <- rep(1, nrow(s$edge))
  weak <- which(sup < 95) + Ntip(s)
  s$edge.length[s$edge[,2] %in% weak] <- 0
  s2 <- di2multi(s, tol=1e-8); s2$node.label <- NULL; s2$edge.length <- NULL
  write.tree(unroot(s2), file.path(O, sprintf("reconciliation/%s.constraint_sp95.nwk", nm)))
  cat(nm, ": tips", Ntip(s), "; species nodes kept (UFBoot>=95)", s2$Nnode, "of", s$Nnode, "\n")
  if (grepl("_all", nm)) {
    f <- intersect(tips, fsm); o <- setdiff(tips, fsm)
    con <- paste0("((", paste(f, collapse=","), "),", paste(o, collapse=","), ");")
    writeLines(con, file.path(O, sprintf("reconciliation/%s.constraint_FSMmono.nwk", nm)))
    cat("  FSM", length(f), "nonFSM", length(o), "\n")
  }
}
