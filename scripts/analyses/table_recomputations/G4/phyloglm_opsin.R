suppressPackageStartupMessages({library(ape);library(phylolm)})
P="/path/to/project/scripts/inputs/"
tr <- read.tree(paste0(P,"myco_bac120.rooted.treefile"))
d <- read.delim(paste0(P,"tsv/S1_genome_metadata.tsv"), stringsAsFactors=FALSE, na.strings=c("NA","","NaN","None"))
d$cl <- ifelse(grepl("^FSM",d$clade),"FSM", ifelse(grepl("^RGM",d$clade),"RGM", ifelse(d$clade=="SGM","SGM",NA)))
d <- d[!is.na(d$cl) & !duplicated(d$tip) & d$tip %in% tr$tip.label,]
dd <- d[d$cl %in% c("FSM","RGM") & !is.na(d$rhodopsin),]; dd$y <- as.integer(dd$rhodopsin=="present")
tt <- drop.tip(tr, setdiff(tr$tip.label, dd$tip)); rownames(dd) <- dd$tip; dd <- dd[tt$tip.label,]
dd$clade <- factor(dd$cl, levels=c("RGM","FSM"))
tab <- table(dd$clade, dd$y); print(tab); ft <- fisher.test(tab); print(ft$p.value); print(ft$estimate)
f <- phyloglm(y ~ clade, data=dd, phy=tt, method="logistic_MPLE", btol=50)
cf <- summary(f)$coefficients; print(cf); cat("alpha", f$alpha, "\n")
