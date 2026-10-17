source("common.R"); suppressPackageStartupMessages(library(PhylogeneticEM))
f <- commandArgs(TRUE)[1]; x <- readRDS(f); r <- x$res; tt <- x$tree; y <- x$Y[1,]
lg <- grepl("_log", f)
out <- data.frame()
for (ms in c("DDSE","Djump","BGHuni","LINselect")) {
  pp <- suppressWarnings(try(params_process(r, method.selection=ms), silent=TRUE)); if (inherits(pp,"try-error")) next
  eds <- pp$shifts$edges; vals <- as.numeric(pp$shifts$values); a <- as.numeric(pp$selection.strength); if (!length(a)) a <- NA
  root <- if (!is.null(pp$root.state$exp.root)) as.numeric(pp$root.state$exp.root) else as.numeric(pp$root.state$value.root)
  # expected tip values under the selected model
  ex <- try(as.numeric(imputed_traits(r, trait=1, save_all=FALSE, where="tips", what="expectations", params=pp)), silent=TRUE)
  for (i in seq_along(eds)) {
    de <- describe_edge(tt, eds[i], y)
    el <- tt$edge.length[eds[i]]
    tips <- if (de$node<=Ntip(tt)) tt$tip.label[de$node] else extract.clade(tt,de$node)$tip.label
    # effective shift in expected tip mean: value*(1-exp(-a*(time from shift to tips))) for OU; value for BM
    dt <- max(node.depth.edgelength(tt)) - node.depth.edgelength(tt)[tt$edge[eds[i],1]]
    eff <- if (is.na(a) || a==0) vals[i] else vals[i]*(1-exp(-a*dt))
    out <- rbind(out, cbind(file=f, method=ms, K=length(eds), alpha=signif(a,4), half_life=signif(log(2)/a,3), root_mean=signif(root,4), de, edge_length=signif(el,3),
      shift_value=signif(vals[i],4), effective_tip_shift=signif(eff,4), effective_ratio=if (lg) signif(exp(eff),3) else NA))
  }
}
print(out[,c("method","K","alpha","n_tips","composition","label","mean_size_Mb_below","edge_length","shift_value","effective_tip_shift","effective_ratio")])
write.table(out, sub("\\.rds$","_shifts.tsv",f), sep="\t", quote=FALSE, row.names=FALSE)
