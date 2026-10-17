#!/usr/bin/env Rscript
# Prevalence-conditioned Mk null (manuscript design, Table S12e: FSM+RGM tree, fitted ARD model, root from the
# stationary distribution, replicates kept if positives are within 10% of observed) run separately per opsin lineage
suppressPackageStartupMessages({library(ape);library(phytools);library(castor)})
# simulation with castor::simulate_mk_model (root from the stationary distribution), in batches of 2,000
O <- "/path/to/project/analyses/extended/B_opsin"
P <- "/path/to/project/scripts"
tr <- read.tree(file.path(P,"inputs/myco_bac120.rooted.treefile"))
d <- read.delim(file.path(P,"inputs/tsv/S1_genome_metadata.tsv"), stringsAsFactors=FALSE, na.strings=c("NA","","NaN","None"))
d$cl <- ifelse(grepl("^FSM",d$clade),"FSM", ifelse(grepl("^RGM",d$clade),"RGM", ifelse(d$clade=="SGM","SGM",NA)))
d$key <- sub("\\.fa$","",d$tip)
gl <- read.delim(file.path(O,"species/genome_lineages.tsv"), stringsAsFactors=FALSE)
d$any_S1 <- as.integer(d$rhodopsin=="present"); for (v in c("any","A","B")) d[[v]] <- as.integer(d$key %in% gl$genome[gl[[v]]==1])
or_of <- function(cl, y){ tab <- table(factor(cl,levels=c("RGM","FSM")), factor(y,levels=c(0,1)))+0.5; (tab[2,2]*tab[1,1])/(tab[2,1]*tab[1,2]) }
dd <- d[d$cl %in% c("FSM","RGM") & d$tip %in% tr$tip.label,]
tt <- drop.tip(tr, setdiff(tr$tip.label, dd$tip)); rownames(dd) <- dd$tip; dd <- dd[tt$tip.label,]
res <- list()
for (v in c("any_S1","any","A","B")) {
  set.seed(20260901)
  yv <- dd[[v]]; y <- setNames(as.character(yv), dd$tip); npos <- sum(yv); obs <- or_of(dd$cl, yv)
  fit <- fitMk(tt, y, model="ARD"); Q <- as.Qmatrix(fit)
  pi1 <- Q["0","1"]/(Q["0","1"]+Q["1","0"])
  lo <- floor(npos*0.9); hi <- ceiling(npos*1.1); acc <- c(); tried <- 0
  Qc <- Q[c("0","1"),c("0","1")]; isF <- dd$cl=="FSM"
  while (length(acc) < 2000 && tried < 500000) {
    sim <- simulate_mk_model(tt, Q=Qc, root_probabilities="stationary", include_nodes=FALSE, Nsimulations=2000, drop_dims=FALSE)$tip_states - 1
    tried <- tried + 2000
    np <- rowSums(sim); ok <- which(np >= lo & np <= hi)
    if (length(ok)) acc <- c(acc, apply(sim[ok,,drop=FALSE], 1, function(s) or_of(dd$cl, s))) }
  acc <- head(acc, 2000)
  p <- (sum(acc>=obs)+1)/(length(acc)+1)
  cat(sprintf("%s: FSM %d/%d, RGM %d/%d; q01 %.3f q10 %.3f; tried %d accepted %d; null OR median %.2f p95 %.2f; observed %.1f; p = %.4f\n",
     v, sum(yv[dd$cl=="FSM"]), sum(dd$cl=="FSM"), sum(yv[dd$cl=="RGM"]), sum(dd$cl=="RGM"), Q["0","1"], Q["1","0"], tried, length(acc), median(acc), quantile(acc,.95), obs, p))
  res[[v]] <- data.frame(trait=v, n=nrow(dd), FSM_pos=sum(yv[dd$cl=="FSM"]), RGM_pos=sum(yv[dd$cl=="RGM"]), q01=Q["0","1"], q10=Q["1","0"],
     tried=tried, accepted=length(acc), null_OR_median=median(acc), null_OR_p95=quantile(acc,.95), observed_OR=obs, p=p)
}
write.table(do.call(rbind,res), file.path(O,"prevalence/prevalence_null_per_lineage.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
