#!/usr/bin/env Rscript
# Heliorhodopsin-rooted tree with SH-aLRT/UFBoot labels re-attached by bipartition (not by node number)
suppressPackageStartupMessages(library(ape))
O <- "/path/to/project/analyses/extended/B_opsin"
u <- read.tree(file.path(O,"tree/opsin_ml.treefile")); r <- read.tree(file.path(O,"tree/opsin_ml.rooted_HeR.nwk"))
lab <- u$tip.label; n <- Ntip(u)
key <- function(s) { s <- sort(s); c2 <- sort(setdiff(lab, s)); if (length(c2) < length(s) || (length(c2)==length(s) && c2[1] < s[1])) s <- c2; paste(s, collapse=",") }
sup <- list(); for (i in (n+1):(n+u$Nnode)) { s <- extract.clade(u, i)$tip.label; if (length(s) < n) sup[[key(s)]] <- u$node.label[i-n] }
r$node.label <- sapply((n+1):(n+r$Nnode), function(i) { s <- extract.clade(r,i)$tip.label; if (length(s) == n) "" else { k <- key(s); if (is.null(sup[[k]])) "" else sup[[k]] } })
write.tree(r, file.path(O,"tree/opsin_ml.rooted_HeR.support.nwk")); cat("labelled", sum(r$node.label != ""), "of", r$Nnode, "nodes\n")
