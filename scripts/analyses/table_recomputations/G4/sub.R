suppressPackageStartupMessages({library(ape)})
P="/path/to/project/scripts/inputs/"
tr <- read.tree(paste0(P,"myco_bac120.rooted.treefile"))
d <- read.delim(paste0(P,"tsv/S1_genome_metadata.tsv"), stringsAsFactors=FALSE, na.strings=c("NA","","NaN","None"))
d$cl <- ifelse(grepl("^FSM",d$clade),"FSM", ifelse(grepl("^RGM",d$clade),"RGM", ifelse(d$clade=="SGM","SGM",NA)))
d <- d[!is.na(d$cl) & !duplicated(d$tip) & d$tip %in% tr$tip.label,]; d$gs <- as.numeric(d$genome_size_Mb)
cat("ingroup n", nrow(d), table(d$cl), "\n")
ing <- drop.tip(tr, setdiff(tr$tip.label, d$tip)); rownames(d) <- d$tip
fsm <- d$tip[d$cl=="FSM"]; cat("FSM monophyletic:", is.monophyletic(ing, fsm), "\n")
ft <- keep.tip(ing, fsm); small <- d[ft$tip.label,"gs"] < 4
cat("FSM small", sum(small), "large", sum(!small), "\n")
# largest all-small clade in a tree given small flag by tip label
largest <- function(tt, sm) {
  names(sm) <- NULL
  pp <- prop.part(tt); best <- if (any(sm)) 1 else 0
  for (i in seq_along(pp)) if (all(sm[pp[[i]]])) best <- max(best, length(pp[[i]]))
  best }
smallv <- setNames(d[ft$tip.label,"gs"]<4, ft$tip.label)
cat("full FSM largest block:", largest(ft, smallv[ft$tip.label]), "\n")
# non-FSM largest block / blocks at 4 Mb
oth <- d$tip[d$cl!="FSM"]
ot <- ing  # blocks of non-FSM small: node all tips non-FSM & small
good <- setNames(d[ing$tip.label,"cl"]!="FSM" & d[ing$tip.label,"gs"]<4, ing$tip.label)
pp <- prop.part(ing); nb <- c(); inb <- rep(FALSE, length(pp))
mx <- 0; for (i in seq_along(pp)) if (all(good[pp[[i]]])) mx <- max(mx, length(pp[[i]]))
cat("largest non-FSM block:", mx, "\n")
set.seed(20261017)
ns <- c(1,2,3,4,5,6,8,10,12,16,20,30,40,60); R <- 2000
cat("n median min max q05 q25 q75 exceed3 P_exceed3\n")
res <- NULL
for (n in ns) { v <- replicate(R, { s <- sample(ft$tip.label, n); if (n==1) 1 else { t2 <- keep.tip(ft, s); largest(t2, smallv[t2$tip.label]) } })
  res <- rbind(res, c(n, median(v), min(v), max(v), quantile(v,.05), quantile(v,.25), quantile(v,.75), sum(v>3), mean(v>3)))
  cat(res[nrow(res),], "\n") }
write.table(res, "rerun/subsample_res.tsv", sep="\t", row.names=FALSE, col.names=c("n","median","min","max","q05","q25","q75","exceed3","p_exceed3"))
# also with ties for 'n genomes needed': P(exceed3) for n=4..20 each
