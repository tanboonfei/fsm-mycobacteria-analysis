suppressPackageStartupMessages({library(ape);library(phylolm)})
tr <- read.tree("/path/to/project/scripts/inputs/myco_bac120.rooted.treefile")
d <- read.delim("fsm.tsv", stringsAsFactors=FALSE, na.strings=c("NA","","NaN","None"))
d <- d[!duplicated(d$tip) & d$tip %in% tr$tip.label,]
cat("n",nrow(d),"\n")
d$opsin <- as.integer(d$rhodopsin=="present"); d$gs <- as.numeric(d$genome_size_Mb); d$comp <- as.numeric(d$completeness_used_pct)
ft <- drop.tip(tr, setdiff(tr$tip.label, d$tip)); rownames(d) <- d$tip; d <- d[ft$tip.label,]
d$ratio <- d$gs/(d$comp/100)
d$ratio_z <- as.numeric(scale(d$ratio)); d$gs_z <- as.numeric(scale(d$gs))
for (f in c("opsin ~ gs","opsin ~ ratio_z","opsin ~ ratio","opsin ~ gs_z")) {
 p <- phyloglm(as.formula(f), data=d, phy=ft, method="logistic_MPLE", btol=30)
 cat(f,"\n"); print(summary(p)$coefficients); cat("alpha",p$alpha,"\n")
}
cat("sd ratio",sd(d$ratio),"sd gs",sd(d$gs),"\n")
