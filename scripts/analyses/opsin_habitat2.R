suppressPackageStartupMessages({library(ape);library(phylolm)})
tr <- read.tree("../Publication/scripts/inputs/myco_bac120.rooted.treefile")
d  <- read.delim("../Publication/scripts/inputs/tsv/S1_genome_metadata.tsv", stringsAsFactors=FALSE, na.strings=c("NA","","NaN","None"))
d <- d[grepl("^FSM",d$clade) & !duplicated(d$tip) & d$tip %in% tr$tip.label,]
d$opsin <- as.integer(d$rhodopsin=="present"); d$gs <- as.numeric(d$genome_size_Mb)
cl <- read.delim("unassigned_ani/clusters.tsv", stringsAsFactors=FALSE)
d$sp <- ifelse(d$species=="novel sp.", cl$cluster[match(d$tip, cl$tip)], d$species)
sink("opsin_habitat2_out.txt", split=TRUE)
f <- d[d$biome %in% c("aquatic","anthropogenic"),]; f$aq <- as.integer(f$biome=="aquatic")
cat("species (ANI-clustered) among anthropogenic:", length(unique(f$sp[f$aq==0])), "; aquatic:", length(unique(f$sp[f$aq==1])), "\n")
print(table(f$sp[f$aq==0], f$opsin[f$aq==0]))
set.seed(1); res <- c(); pv <- c()
for (rep in 1:200) {
  g <- do.call(rbind, lapply(split(f, f$sp), function(z) z[sample(nrow(z),1),]))
  gt <- drop.tip(tr, setdiff(tr$tip.label, g$tip)); rownames(g) <- g$tip; g <- g[gt$tip.label,]
  pv <- c(pv, fisher.test(table(g$aq, g$opsin))$p.value)
  if (rep==1) { cat("one genome per species, n =", nrow(g), "\n"); print(table(aquatic=g$aq, opsin=g$opsin))
    p <- suppressWarnings(phyloglm(opsin ~ aq, data=g, phy=gt, method="logistic_MPLE", btol=30)); print(summary(p)$coefficients); cat("alpha", p$alpha, "\n") }
}
cat("Fisher p across 200 random one-per-species draws: median", signif(median(pv),3), " 95th pct", signif(quantile(pv,.95),3), "\n")
sink()
