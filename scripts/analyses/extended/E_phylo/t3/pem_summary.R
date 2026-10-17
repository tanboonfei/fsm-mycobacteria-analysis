source("common.R"); suppressPackageStartupMessages(library(PhylogeneticEM))
fs <- commandArgs(TRUE); out <- data.frame(); path <- data.frame()
for (f in fs) {
  x <- readRDS(f); r <- x$res; tt <- x$tree
  fe <- which(tt$edge[,2]==getMRCA(tt, intersect(fsm_tips, tt$tip.label)))
  lp <- intersect(lep_tips, tt$tip.label); m <- getMRCA(tt, lp); le <- c()
  while (TRUE) { nd <- if (m<=Ntip(tt)) 1 else length(extract.clade(tt,m)$tip.label); if (nd>5) break; le <- c(le, which(tt$edge[,2]==m)); m <- tt$edge[tt$edge[,2]==m,1]; if (m==Ntip(tt)+1) break }
  for (ms in c("DDSE","Djump","BGHuni")) {
    pp <- suppressWarnings(try(params_process(r, method.selection=ms), silent=TRUE)); if (inherits(pp,"try-error")) next
    e <- pp$shifts$edges
    out <- rbind(out, data.frame(run=sub("\\.rds$","",f), tips=Ntip(tt), process=r$process, method=ms, K_selected=length(e), alpha=signif(as.numeric(pp$selection.strength),3),
      FSM_stem_shift=fe %in% e, leprosy_clade_shift=any(le %in% e), n_single_tip_shifts=sum(tt$edge[e,2] <= Ntip(tt))))
  }
  # best configuration at each K (alpha_max), does it contain the FSM stem?
  pe <- r$alpha_max$params_estim
  for (k in seq_along(pe)) { e <- pe[[k]]$shifts$edges; ll <- attr(pe[[k]], "log_likelihood")
    path <- rbind(path, data.frame(run=sub("\\.rds$","",f), K=k-1, logLik=if (is.null(ll)) NA else round(as.numeric(ll[1]),2), FSM_stem=fe %in% e, leprosy=any(le %in% e))) }
}
print(out); write.table(out, "pem_selection_summary.tsv", sep="\t", quote=FALSE, row.names=FALSE)
w <- tapply(seq_len(nrow(path)), path$run, function(i) { p <- path[i,]; c(first_K_with_FSM=if (any(p$FSM_stem)) min(p$K[p$FSM_stem]) else NA, first_K_with_leprosy=if (any(p$leprosy)) min(p$K[p$leprosy]) else NA) })
print(w); write.table(path, "pem_path_byK.tsv", sep="\t", quote=FALSE, row.names=FALSE)
