## Stochastic-mapping check (phytools make.simmap, 100 maps per family, Q fixed at the model value, pi = stationary)
## of the analytic joint P(1->0) and expected loss events on the FSM stem edge, for the 53 clade-wide-undetected families.
## make.simmap samples edge histories by rejection, which can stall on very short edges whose sampled end states differ,
## so each family gets a 180-s limit and results are appended per family.  usage: Rscript simmap_check.R <panel> <model>
suppressPackageStartupMessages({library(ape); library(phytools)})
args <- commandArgs(TRUE); pn <- args[1]; mod <- args[2]
DST <- "/path/to/project/analyses/relax_stem/stem/out/"
d <- read.delim(gzfile(paste0(DST, "joint_stem_per_family.tsv.gz")))
set.seed(20261016)
P <- readRDS(paste0(DST, "panel_", pn, ".rds")); tt <- P$tt; k_stem <- which(tt$edge[, 2] == P$crown)
e <- d[d$panel == pn & d$model == mod & d$fate == "lost_cladewide", ]
of <- paste0(DST, "simmap_check_53_", pn, "_", gsub("[^A-Za-z]", "", mod), ".tsv")
cat("panel\tmodel\tfamily\tgain\tjoint_p10_analytic\tjoint_p10_simmap\tE_loss_analytic\tE_loss_simmap\tsecs\n", file = of)
for (i in seq_len(nrow(e))) {
  fam <- e$family[i]; g <- e$gain[i]; l <- e$loss[i]; t0 <- Sys.time(); r <- c(NA, NA)
  if (g <= 50) {
    x <- setNames(as.character(P$X[, fam]), rownames(P$X))
    Q <- matrix(c(-g, g, l, -l), 2, 2, byrow = TRUE, dimnames = list(c("0","1"), c("0","1")))
    sm <- tryCatch({ setTimeLimit(elapsed = 180, transient = TRUE)
      make.simmap(tt, x, Q = Q, pi = c(l, g) / (g + l), nsim = 100, message = FALSE) }, error = function(err) NULL)
    setTimeLimit(elapsed = Inf)
    if (!is.null(sm)) r <- rowMeans(sapply(sm, function(s) { st <- names(s$maps[[k_stem]]); c(st[1] == "1" & tail(st, 1) == "0", sum(head(st, -1) == "1" & tail(st, -1) == "0")) }))
  }
  cat(pn, mod, fam, g, e$joint_p10[i], r[1], e$E_loss_events_stem[i], r[2], round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), sep = "\t", file = of, append = TRUE); cat("\n", file = of, append = TRUE)
}
