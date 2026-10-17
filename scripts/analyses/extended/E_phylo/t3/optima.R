source("common.R"); suppressPackageStartupMessages(library(l1ou))
f <- commandArgs(TRUE)[1]; x <- readRDS(f); fit <- x$fit; tt <- x$tree
cat(names(fit), "\n")
op <- fit$optima; if (is.matrix(op)) op <- op[,1]; names(op) <- tt$tip.label
lg <- grepl("_log_", f); tr <- function(v) if (lg) exp(v) else v
grp <- list(FSM=intersect(fsm_tips,tt$tip.label), sister42=intersect(s42,tt$tip.label), leprosy=intersect(lep_tips,tt$tip.label))
for (g in names(grp)) cat(g, "optimum(s):", paste(unique(round(tr(op[grp[[g]]]),3)), collapse=","), " modal:", round(tr(as.numeric(names(which.max(table(round(op[grp[[g]]],5)))))),3), "\n")
cat("intercept (root optimum):", round(tr(fit$intercept),3), "\n")
