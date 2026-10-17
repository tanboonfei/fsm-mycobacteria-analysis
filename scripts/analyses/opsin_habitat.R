suppressPackageStartupMessages({library(ape);library(phylolm)})
tr <- read.tree("../Publication/scripts/inputs/myco_bac120.rooted.treefile")
d  <- read.delim("../Publication/scripts/inputs/tsv/S1_genome_metadata.tsv", stringsAsFactors=FALSE, na.strings=c("NA","","NaN","None"))
d <- d[grepl("^FSM",d$clade) & !duplicated(d$tip) & d$tip %in% tr$tip.label,]
d$opsin <- as.integer(d$rhodopsin=="present"); d$gs <- as.numeric(d$genome_size_Mb)
sink("opsin_habitat_out.txt", split=TRUE)
cat("FSM opsin by biome:\n"); print(table(d$biome, d$opsin, useNA="ifany"))
f <- d[d$biome %in% c("aquatic","anthropogenic"),]; f$aq <- as.integer(f$biome=="aquatic")
ft <- drop.tip(tr, setdiff(tr$tip.label, f$tip)); rownames(f) <- f$tip; f <- f[ft$tip.label,]
ft <- compute.brlen(ft, method=function(n) ft$edge.length) # keep lengths
ft$edge.length[ft$edge.length<=0] <- 1e-6
ft2 <- read.tree(text=write.tree(drop.tip(tr, setdiff(tr$tip.label, f$tip))))
tb <- table(f$aq, f$opsin); print(tb); print(fisher.test(tb))
cat("n species among anthropogenic:", length(unique(f$species[f$aq==0])), "; aquatic:", length(unique(f$species[f$aq==1])), "\n")
print(f[f$aq==0, c("tip","species","gs","opsin")])
p1 <- phyloglm(opsin ~ aq, data=f, phy=ft2, method="logistic_MPLE", btol=30); cat("\nphyloglm opsin ~ aquatic\n"); print(summary(p1)$coefficients); cat("alpha", p1$alpha, "\n")
p2 <- phyloglm(opsin ~ aq + gs, data=f, phy=ft2, method="logistic_MPLE", btol=30); cat("\nphyloglm opsin ~ aquatic + genome size\n"); print(summary(p2)$coefficients); cat("alpha", p2$alpha, "\n")
# one genome per species (collapse; keep unnamed as separate)
f$sp <- ifelse(f$species=="novel sp.", f$tip, f$species)
g <- f[!duplicated(f$sp),]; gt <- drop.tip(ft2, setdiff(ft2$tip.label, g$tip)); rownames(g) <- g$tip; g <- g[gt$tip.label,]
cat("\nOne genome per species: n =", nrow(g)); print(table(g$aq, g$opsin)); print(fisher.test(table(g$aq,g$opsin))$p.value)
p3 <- tryCatch(phyloglm(opsin ~ aq, data=g, phy=gt, method="logistic_MPLE", btol=30), error=function(e) NULL); if(!is.null(p3)) print(summary(p3)$coefficients)
sink()
