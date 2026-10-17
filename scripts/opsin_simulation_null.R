#!/usr/bin/env Rscript
# Phylogenetic simulation null for the opsin association, matching the design
# already used for aquatic biome and low salinity (revision_checks/m6c_simnull.R).
suppressPackageStartupMessages({library(ape); library(phytools)})
set.seed(20260901)
tr <- read.tree("inputs/myco_bac120.rooted.treefile")
d  <- read.delim("inputs/tsv/S1_genome_metadata.tsv",
                 stringsAsFactors = FALSE, na.strings = c("NA","","NaN","None"))
d$cl <- ifelse(grepl("^FSM", d$clade), "FSM",
        ifelse(grepl("^RGM", d$clade), "RGM",
        ifelse(d$clade == "SGM", "SGM", NA)))
d$opsin <- ifelse(is.na(d$rhodopsin), NA, as.integer(d$rhodopsin == "present"))
or_of <- function(clade, y) {
  tab <- table(factor(clade, levels=c("RGM","FSM")), factor(y, levels=c(0,1))) + 0.5
  (tab[2,2]*tab[1,1])/(tab[2,1]*tab[1,2])
}
SIMS <- 2000
dd <- d[d$cl %in% c("FSM","RGM") & !is.na(d$opsin), ]
dd <- dd[dd$tip %in% tr$tip.label, ]
tt <- drop.tip(tr, setdiff(tr$tip.label, dd$tip))
rownames(dd) <- dd$tip; dd <- dd[tt$tip.label, ]
y <- dd$opsin; names(y) <- dd$tip
fit <- fitMk(tt, setNames(as.character(y), names(y)), model="ER")
rate <- fit$rates[1]
obs <- or_of(dd$cl, y)
cat(sprintf("n = %d (FSM %d, RGM %d); observed OR = %.2f; ER rate = %.4f\n",
            length(y), sum(dd$cl=="FSM"), sum(dd$cl=="RGM"), obs, rate))
ors <- numeric(0)
for (i in seq_len(SIMS)) {
  s <- tryCatch(rTraitDisc(tt, model="ER", k=2, rate=rate, states=c("0","1"),
                           root.value=sample(1:2,1)), error=function(e) NULL)
  if (is.null(s)) next
  sv <- as.integer(as.character(s[tt$tip.label]))
  ors <- c(ors, or_of(dd$cl, sv))
}
p <- (sum(ors >= obs) + 1)/(length(ors) + 1)
cat(sprintf("usable sims = %d; sim OR median = %.2f, p95 = %.2f; simulation-null p = %.4f\n",
            length(ors), median(ors), quantile(ors,0.95), p))
dir.create("run_output", showWarnings=FALSE)
write.table(data.frame(test="Opsin presence", reference="RGM", n=length(y),
                       observed_OR=obs, sim_OR_median=median(ors),
                       sim_OR_p95=as.numeric(quantile(ors,0.95)),
                       simulation_null_p=p, n_sims=length(ors)),
            "run_output/S3b_opsin_simnull.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)
