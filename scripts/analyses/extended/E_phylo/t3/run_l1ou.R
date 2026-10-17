source("common.R"); suppressPackageStartupMessages(library(l1ou))
args <- commandArgs(TRUE); set <- args[1]; scale <- args[2]; crit <- args[3]; treef <- args[4]; tree_type <- sub("\\.nwk$","",sub("myco_bac120\\.","",basename(treef)))
groups <- if (set=="824") c("FSM","RGM","SGM") else c("FSM","RGM")
derep <- length(args)>4 && args[5]=="derep"; if (derep) tree_type <- paste0(tree_type,"_derep")
tu <- read.tree(treef); keep <- if (derep) intersect(derep_tips(groups), tu$tip.label) else d$tip[d$cl %in% groups & d$tip %in% tu$tip.label]
tt <- drop.tip(tu, setdiff(tu$tip.label, keep)); tt$node.label <- NULL; tt <- reorder(tt,"postorder")
y <- d[tt$tip.label,"genome_size_Mb"]; if (scale=="log") y <- log(y)
names(y) <- tt$tip.label
ad <- adjust_data(tt, y); tt <- ad$tree; Y <- ad$Y
t0 <- Sys.time()
fit <- estimate_shift_configuration(tt, Y, criterion=crit, max.nShifts=25, nCores=1)
cat("time", format(Sys.time()-t0), "\n"); print(fit)
sh <- fit$shift.configuration
res <- do.call(rbind, lapply(seq_along(sh), function(i) cbind(describe_edge(tt, sh[i], Y), shift_value=round(fit$shift.values[i],4))))
res <- cbind(set=set, scale=scale, criterion=crit, tree=tree_type, alpha=signif(fit$alpha,4), sigma2=signif(fit$sigma2,4), intercept=signif(fit$intercept,4), nShifts=length(sh), score=round(fit$score,2), res)
print(res[,c("n_tips","composition","label","mean_size_Mb_below","shift_value","species_examples")])
write.table(res, sprintf("l1ou_%s_%s_%s_%s.tsv", set, scale, crit, tree_type), sep="\t", quote=FALSE, row.names=FALSE)
saveRDS(list(fit=fit, tree=tt, Y=Y), sprintf("l1ou_%s_%s_%s_%s.rds", set, scale, crit, tree_type))
