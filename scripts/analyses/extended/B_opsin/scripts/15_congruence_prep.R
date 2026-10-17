#!/usr/bin/env Rscript
# Gene-tree vs species-tree congruence per lineage: RF with label-permutation null, strongly supported conflicts,
# and inputs for the AU test and GeneRax.
suppressPackageStartupMessages({library(ape); library(phangorn)})
set.seed(20261014)
O <- "/path/to/project/analyses/extended/B_opsin"
P <- "/path/to/project/scripts"
u <- read.tree(file.path(O,"tree/opsin_ml.treefile"))
sp <- read.tree(file.path(P,"inputs/myco_bac120.rooted.treefile")); sp$node.label <- NULL
L <- read.delim(file.path(O,"screen/opsin_loci_lineages.tsv"), stringsAsFactors=FALSE)
L <- L[L$in_tree %in% c("True","TRUE",TRUE), ]
L$sptip <- ifelse(grepl("^spire_mag", L$genome), paste0(L$genome, ".fa"), NA)
s1 <- read.delim(file.path(P,"inputs/tsv/S1_genome_metadata.tsv"), stringsAsFactors=FALSE)
s1$key <- sub("\\.fa$","",s1$tip); L$sptip <- s1$tip[match(L$genome, s1$key)]
stopifnot(all(L$sptip %in% sp$tip.label))
aln <- read.FASTA(file.path(O,"tree/opsins_all.trim.aln"), type="AA")
res <- list()
sets <- list(A_all = L[L$lineage=="A",], A_FSM = L[L$lineage=="A" & L$clade=="FSM",], B_all = L[L$lineage=="B",], B_FSM = L[L$lineage=="B" & L$clade=="FSM",])
for (nm in names(sets)) {
  x <- sets[[nm]]; x <- x[order(-x$aa_len), ]; ncopy <- sum(duplicated(x$genome)); x <- x[!duplicated(x$genome), ]
  g <- keep.tip(u, x$locus); g$tip.label <- x$sptip[match(g$tip.label, x$locus)]
  s <- unroot(keep.tip(sp, x$sptip)); g <- unroot(g)
  rf <- RF.dist(g, s, normalize=TRUE)
  perm <- replicate(1000, { s2 <- s; s2$tip.label <- sample(s2$tip.label); RF.dist(g, s2, normalize=TRUE) })
  # strongly supported gene-tree splits that conflict with the species tree
  gs <- keep.tip(u, x$locus); n <- Ntip(gs); sup <- gs$node.label
  ufb <- suppressWarnings(as.numeric(sub(".*/","",sup))); shl <- suppressWarnings(as.numeric(sub("/.*","",sup)))
  gs$tip.label <- x$sptip[match(gs$tip.label, x$locus)]
  sps <- as.splits(s); comp <- 0; conf <- 0; strong <- 0; strong_conf <- 0
  for (i in (n+2):(n+gs$Nnode)) { tips <- extract.clade(gs, i)$tip.label; if (length(tips) < 2 || length(tips) > n-2) next
    ok <- is.monophyletic(s, tips) || is.monophyletic(s, setdiff(gs$tip.label, tips))
    st <- !is.na(ufb[i-n]) && ufb[i-n] >= 95 && !is.na(shl[i-n]) && shl[i-n] >= 80
    if (ok) comp <- comp+1 else conf <- conf+1; if (st) { strong <- strong+1; if (!ok) strong_conf <- strong_conf+1 } }
  res[[nm]] <- data.frame(set=nm, genomes=nrow(x), extra_copies_dropped=ncopy, nRF=round(rf,3), perm_nRF_median=round(median(perm),3),
     perm_nRF_min=round(min(perm),3), p_perm=(sum(perm<=rf)+1)/1001, gene_splits_compatible=comp, gene_splits_conflicting=conf,
     strong_splits=strong, strong_splits_conflicting=strong_conf)
  # AU-test inputs: sub-alignment with species-tip names, species topology
  sub <- aln[x$locus]; names(sub) <- x$sptip; write.FASTA(sub, file.path(O, sprintf("reconciliation/%s.aln", nm)))
  write.tree(s, file.path(O, sprintf("reconciliation/%s.speciestopo.nwk", nm)))
  write.tree(g, file.path(O, sprintf("reconciliation/%s.genetree_pruned.nwk", nm)))
}
R <- do.call(rbind, res); print(R)
write.table(R, file.path(O,"reconciliation/rf_congruence.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
# GeneRax inputs: all mycobacterial copies per lineage
spf <- sp; write.tree(spf, file.path(O,"reconciliation/species_862.nwk"))
fam <- c("[FAMILIES]")
for (lin in c("A","B")) {
  x <- L[L$lineage==lin,]
  g <- unroot(keep.tip(u, x$locus)); g$node.label <- NULL; write.tree(g, file.path(O, sprintf("reconciliation/gr_%s.start.nwk", lin)))
  sub <- aln[x$locus]; write.FASTA(sub, file.path(O, sprintf("reconciliation/gr_%s.aln", lin)))
  writeLines(c("[ALL]", "", paste0(x$locus, ":", x$sptip)), file.path(O, sprintf("reconciliation/gr_%s.map", lin)))
  # GeneRax link format: one line per gene "gene:species"
  writeLines(paste(x$locus, x$sptip, sep="\t"), file.path(O, sprintf("reconciliation/gr_%s.map", lin)))
  fam <- c(fam, paste0("- lineage", lin), paste0("starting_gene_tree = ", O, "/reconciliation/gr_", lin, ".start.nwk"),
           paste0("alignment = ", O, "/reconciliation/gr_", lin, ".aln"), paste0("mapping = ", O, "/reconciliation/gr_", lin, ".map"), "subst_model = LG+G4")
}
writeLines(fam, file.path(O,"reconciliation/families.txt"))
