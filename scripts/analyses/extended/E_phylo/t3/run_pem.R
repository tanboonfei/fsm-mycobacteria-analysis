# PhylogeneticEM shift detection. args: set scale treefile ncores process(scOU|BM) derep(0|1) me(0|1) tag
source("common.R"); suppressPackageStartupMessages({library(PhylogeneticEM); library(phylolm)})
a <- commandArgs(TRUE); set <- a[1]; scale <- a[2]; treef <- a[3]; ncores <- as.integer(a[4]); proc <- a[5]; derep <- a[6]=="1"; me <- a[7]=="1"; tag <- a[8]
groups <- if (set=="824") c("FSM","RGM","SGM") else c("FSM","RGM")
keep <- if (derep) derep_tips(groups) else d$tip[d$cl %in% groups]
tu <- read.tree(treef); tt <- drop.tip(tu, setdiff(tu$tip.label, keep)); tt$node.label <- NULL
tt$edge.length <- tt$edge.length / max(node.depth.edgelength(tt))
y <- d[tt$tip.label,"genome_size_Mb"]; if (scale=="log") y <- log(y); names(y) <- tt$tip.label
if (me) { # tip measurement-error variance relative to BM rate, from an intercept-only BM fit with measurement error
  fm <- phylolm(y ~ 1, data=data.frame(y=y, row.names=names(y)), phy=tt, model="BM", measurement_error=TRUE)
  cc <- fm$sigma2_error / fm$sigma2; cat("ME fit: sigma2", fm$sigma2, " sigma2_error", fm$sigma2_error, " added to terminals (tree heights):", cc, "\n")
  term <- tt$edge[,2] <= Ntip(tt); tt$edge.length[term] <- tt$edge.length[term] + cc
  tt$edge.length <- tt$edge.length / max(node.depth.edgelength(tt))
}
if (proc=="scOU") stopifnot(is.ultrametric(tt))
cat("tips", Ntip(tt), " FSM", sum(tt$tip.label %in% fsm_tips), "\n")
Y <- matrix(y, nrow=1, dimnames=list("size", tt$tip.label))
set.seed(7); t0 <- Sys.time()
res <- PhyloEM(phylo=tt, Y_data=Y, process=proc, random.root=(proc=="scOU"), stationary.root=(proc=="scOU"),
               K_max=20, nbr_alpha=8, parallel_alpha=(ncores>1), Ncores=ncores, progress.bar=FALSE)
cat("time", format(Sys.time()-t0), "\n")
saveRDS(list(res=res, tree=tt, Y=Y), sprintf("pem_%s_%s%s.rds", set, scale, tag))
