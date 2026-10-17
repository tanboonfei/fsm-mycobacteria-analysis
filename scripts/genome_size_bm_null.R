suppressPackageStartupMessages({library(ape);library(phytools);library(geiger)})
set.seed(42)
tr <- read.tree("inputs/myco_bac120.rooted.treefile")
d  <- read.delim("inputs/tsv/S1_genome_metadata.tsv", stringsAsFactors=FALSE)
d$cl <- ifelse(grepl("^FSM",d$clade),"FSM", ifelse(grepl("^RGM",d$clade),"RGM",
        ifelse(d$clade=="SGM","SGM",NA)))
run <- function(groups, tag) {
  dd <- d[d$cl %in% groups & !is.na(d$genome_size_Mb),]
  dd <- dd[dd$tip %in% tr$tip.label,]
  tt <- drop.tip(tr, setdiff(tr$tip.label, dd$tip))
  rownames(dd) <- dd$tip; dd <- dd[tt$tip.label,]
  y <- setNames(as.numeric(dd$genome_size_Mb), rownames(dd))
  fsm <- rownames(dd)[dd$cl=="FSM"]; oth <- setdiff(names(y), fsm)
  obs <- mean(y[fsm]) - mean(y[oth])
  fit <- fitContinuous(tt, y, model="BM")
  sig2 <- fit$opt$sigsq
  sims <- fastBM(tt, sig2=sig2, nsim=2000, a=mean(y))
  nul <- apply(sims, 2, function(z) mean(z[fsm]) - mean(z[oth]))
  p <- (sum(abs(nul) >= abs(obs)) + 1)/(length(nul)+1)
  cat(sprintf("%-22s n=%3d obs %+6.2f Mb | BM sigsq %.4f | null sd %.2f | p = %.3f\n",
      tag, length(y), obs, sig2, sd(nul), p))
}
run(c("FSM","RGM"), "FSM vs RGM")
run(c("FSM","RGM","SGM"), "FSM vs full ingroup")

# Note. analysis_phylo/bmsim.R fits an OU model whose alpha rails to the upper
# bound (half-life 0.014), collapsing the null variance to zero and returning
# p = 1/2001 for any observed shift. That fit is reported in Table S23C and
# marked not reportable. The Brownian null above is the one the manuscript
# quotes: null sd 4.78 Mb, p = 0.60 for FSM against RGM.
