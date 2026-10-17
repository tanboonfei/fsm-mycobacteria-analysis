suppressPackageStartupMessages({library(ape); library(phytools)})
tr <- read.tree("myco_bac120.rooted.treefile"); tr$node.label <- NULL
set.seed(1); t0 <- Sys.time()
ch <- chronos(tr, lambda=1, model="correlated", quiet=TRUE); class(ch) <- "phylo"
cat("chronos time", format(Sys.time()-t0), " ultrametric", is.ultrametric(ch), " PHIIC", attr(ch,"PHIIC")$PHIIC, "\n")
write.tree(ch, "myco_bac120.chronos_l1.nwk")
ex <- force.ultrametric(tr, method="extend"); write.tree(ex, "myco_bac120.extend.nwk")
# correlation of node depths with root-to-tip distances
cat("cor(phylogram RTT, ingroup)", "\n")
