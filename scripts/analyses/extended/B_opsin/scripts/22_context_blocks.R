#!/usr/bin/env Rscript
# Maximal context-pure clades within each lineage subtree (loci without context info are ignored)
suppressPackageStartupMessages(library(ape))
O <- "/path/to/project/analyses/extended/B_opsin"
cx <- read.delim(file.path(O,"synteny/opsin_context.tsv"), stringsAsFactors=FALSE); ctx <- setNames(cx$context, cx$locus)
A <- read.delim(file.path(O,"tree/tip_annotation.tsv"), stringsAsFactors=FALSE); grp <- setNames(A$group, A$tip)
for (lin in c("A","B")) {
  st <- read.tree(file.path(O, sprintf("tree/lineage%s_subtree.nwk", lin)))
  st <- keep.tip(st, intersect(st$tip.label, names(ctx)[ctx != "other"]))
  lab <- paste0(ctx[st$tip.label], "_", grp[st$tip.label])
  n <- Ntip(st); nodes <- (n+1):(n+st$Nnode)
  pure <- nodes[sapply(nodes, function(i) length(unique(lab[match(extract.clade(st,i)$tip.label, st$tip.label)]))==1)]
  maxp <- pure[sapply(pure, function(i) !(st$edge[st$edge[,2]==i,1] %in% pure))]
  blocks <- sapply(maxp, function(i) { t <- extract.clade(st,i)$tip.label; paste0(lab[match(t[1], st$tip.label)], "(", length(t), ")") })
  covered <- unlist(lapply(maxp, function(i) extract.clade(st,i)$tip.label))
  sing <- table(lab[!(st$tip.label %in% covered)])
  cat("\nlineage", lin, ": context x clade blocks:", paste(sort(blocks), collapse=" "), "\n  singletons:", paste(names(sing), sing, collapse=" "), "\n")
  # tip order (ladderized) as a compact string
  st2 <- ladderize(st); o <- st2$edge[st2$edge[,2] <= n, 2]
  rl <- rle(lab[match(st2$tip.label[o], st$tip.label)]); cat("  ladderized order:", paste0(rl$values, "x", rl$lengths, collapse=" | "), "\n")
}
