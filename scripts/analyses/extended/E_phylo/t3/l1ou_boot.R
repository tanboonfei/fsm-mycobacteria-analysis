source("common.R"); suppressPackageStartupMessages(library(l1ou))
a <- commandArgs(TRUE); f <- a[1]; n <- as.integer(a[2]); nc <- as.integer(a[3])
x <- readRDS(f); set.seed(11); t0 <- Sys.time()
bs <- l1ou_bootstrap_support(x$fit, nItrs=n, multicore=(nc>1), nCores=nc)
cat("time", format(Sys.time()-t0), "\n")
tt <- x$tree; sup <- bs$detection.rate
fsm_e <- which(tt$edge[,2]==getMRCA(tt, intersect(fsm_tips, tt$tip.label)))
lep <- intersect(lep_tips, tt$tip.label); lep_e <- if (length(lep)) which(tt$edge[,2]==getMRCA(tt, lep)) else NA
cat("FSM stem edge", fsm_e, "support", sup[fsm_e], "\n")
cat("leprae+lepromatosis MRCA edge", lep_e, "support", sup[lep_e], "\n")
sh <- x$fit$shift.configuration
out <- do.call(rbind, lapply(sh, function(e) cbind(describe_edge(tt, e, x$Y[,1]), bootstrap_detection=sup[e])))
# also edges with support >= 0.3 not in the configuration
oth <- setdiff(which(sup >= 0.3), sh); if (length(oth)) out <- rbind(out, do.call(rbind, lapply(oth, function(e) cbind(describe_edge(tt, e, x$Y[,1]), bootstrap_detection=sup[e]))))
print(out[,c("n_tips","composition","label","mean_size_Mb_below","bootstrap_detection")])
write.table(out, sub("\\.rds$", sprintf("_boot%d.tsv", n), f), sep="\t", quote=FALSE, row.names=FALSE)
