#!/usr/bin/env Rscript
# Per-lineage mapping of opsin presence onto the bac120 species tree (824-tip ingroup):
# dispersion, Sankoff parsimony (free and node-constrained), Mk (ER/ARD) with node-conditioned likelihoods,
# stochastic mapping, and per-lineage phyloglm (manuscript design).
suppressPackageStartupMessages({library(ape); library(phytools); library(castor); library(phylolm)})
set.seed(20261014)
O <- "/path/to/project/analyses/extended/B_opsin"
P <- "/path/to/project/scripts"
tr <- read.tree(file.path(P, "inputs/myco_bac120.rooted.treefile"))
d  <- read.delim(file.path(P, "inputs/tsv/S1_genome_metadata.tsv"), stringsAsFactors=FALSE, na.strings=c("NA","","NaN","None"))
d$cl <- ifelse(grepl("^FSM", d$clade), "FSM", ifelse(grepl("^RGM", d$clade), "RGM", ifelse(d$clade=="SGM","SGM",NA)))
d$key <- sub("\\.fa$","",d$tip)
gl <- read.delim(file.path(O, "species/genome_lineages.tsv"), stringsAsFactors=FALSE)
traits <- setdiff(names(gl), "genome")
for (v in traits) d[[v]] <- as.integer(d$key %in% gl$genome[gl[[v]] == 1])
d$any_S1 <- as.integer(d$rhodopsin == "present")
d <- d[!is.na(d$cl) & d$tip %in% tr$tip.label, ]
tt <- drop.tip(tr, setdiff(tr$tip.label, d$tip)); tt$node.label <- NULL
rownames(d) <- d$tip; d <- d[tt$tip.label, ]
N <- Ntip(tt)
fsm <- d$tip[d$cl=="FSM"]; nF <- getMRCA(tt, fsm)
stopifnot(length(extract.clade(tt, nF)$tip.label) == length(fsm))
nFS <- tt$edge[tt$edge[,2]==nF,1]              # FSM + sister MRCA (parent of FSM crown)
nR  <- N + 1                                    # ingroup root
cat(sprintf("ingroup %d tips; FSM %d; FSM+sister clade %d tips\n", N, length(fsm), length(extract.clade(tt,nFS)$tip.label)))
desc <- function(n) if (n <= N) n else unlist(lapply(tt$edge[tt$edge[,1]==n,2], desc))
# ---------- Sankoff parsimony with optional node constraints ----------
tt <- reorder(tt, "cladewise"); ch <- split(tt$edge[,2], tt$edge[,1])
nodes_up <- rev(unique(tt$edge[,1]))   # every node after all its descendants
sankoff <- function(y, gain=1, loss=1, fix=list()) {
  C <- matrix(c(0, gain, loss, 0), 2, byrow=TRUE)  # C[from,to], states 0/1
  S <- matrix(Inf, N + tt$Nnode, 2); S[cbind(1:N, y+1)] <- 0
  nodes <- nodes_up
  for (n in nodes) {
    k <- ch[[as.character(n)]]
    s <- c(0,0)
    for (c in k) s <- s + sapply(1:2, function(a) min(C[a,] + S[c,]))
    if (!is.null(fix[[as.character(n)]])) s[-(fix[[as.character(n)]]+1)] <- Inf
    S[n,] <- s
  }
  min(S[nR,])
}
# count gains/losses on an ACCTRAN-free MPR via castor (Fitch, unordered)
mpr_counts <- function(y, wl=1) {
  tm <- matrix(c(0,1,wl,0),2,byrow=TRUE)
  a <- asr_max_parsimony(tt, y+1, Nstates=2, transition_costs=tm)
  st <- rbind(matrix(0,N,2), a$ancestral_likelihoods); st[cbind(1:N,y+1)] <- 1
  s <- max.col(st, ties.method="first")
  pa <- s[tt$edge[,1]]; cc <- s[tt$edge[,2]]
  list(gains=sum(pa==1 & cc==2), losses=sum(pa==2 & cc==1), states=s)
}
# joint posterior sampling of node states (Felsenstein pruning with scaling; FitzJohn root prior)
pre_edges <- NULL
joint_sample_counts <- function(y, Q, nsamp) {
  nn <- N + tt$Nnode; Lk <- matrix(1, nn, 2); Lk[1:N,] <- 0; Lk[cbind(1:N, y+1)] <- 1
  Pm <- lapply(seq_len(nrow(tt$edge)), function(e) expm::expm(Q * tt$edge.length[e]))
  eidx <- split(seq_len(nrow(tt$edge)), tt$edge[,1])
  for (n in nodes_up) { v <- c(1,1); for (e in eidx[[as.character(n)]]) v <- v * as.vector(Pm[[e]] %*% Lk[tt$edge[e,2],]); Lk[n,] <- v/sum(v) }
  pr <- Lk[nR,]/sum(Lk[nR,])
  inF <- tt$edge[,2] %in% getDescendants(tt, nF); stem <- which(tt$edge[,2]==nF)
  ord <- seq_len(nrow(tt$edge))   # cladewise order: parent edge before child edges
  out <- matrix(0, nsamp, 6, dimnames=list(NULL, c("gains","losses","gains_inFSM","losses_inFSM","gain_on_stem","gains_outFSM")))
  for (k in seq_len(nsamp)) {
    st <- integer(nn); st[nR] <- sample(1:2, 1, prob=pr)
    for (e in ord) { c <- tt$edge[e,2]; if (c <= N) { st[c] <- y[c]+1; next }
      w <- Pm[[e]][st[tt$edge[e,1]],] * Lk[c,]; st[c] <- sample(1:2, 1, prob=w) }
    pa <- st[tt$edge[,1]]; ch2 <- st[tt$edge[,2]]
    g <- pa==1 & ch2==2; l <- pa==2 & ch2==1
    out[k,] <- c(sum(g), sum(l), sum(g[inF]), sum(l[inF]), as.integer(g[stem]), sum(g[!inF & seq_along(g)!=stem]))
  }
  out
}
# ---------- Mk helpers ----------
fit_best <- function(y) {
  x <- setNames(as.character(y), d$tip)
  er <- fitMk(tt, x, model="ER", pi="fitzjohn"); ard <- fitMk(tt, x, model="ARD", pi="fitzjohn")
  aic <- c(ER=AIC(er), ARD=AIC(ard))
  list(ER=er, ARD=ard, aic=aic, best=if (aic["ARD"] + 0 < aic["ER"]) ard else er, bestname=names(which.min(aic)))
}
# likelihood conditioned on node states: attach zero-length pendant tips carrying the state
cond_lik <- function(y, Q, fixnodes) {
  t2 <- tt; x <- setNames(as.character(y), d$tip)
  # attach in decreasing node order so node numbers of later targets are tracked by MRCA labels
  tgt <- lapply(names(fixnodes), function(n) { n <- as.integer(n); if (n==nR) d$tip else extract.clade(tt, n)$tip.label })
  for (i in seq_along(tgt)) {
    nn <- if (length(tgt[[i]]) == Ntip(t2) - (i-1) || identical(sort(tgt[[i]]), sort(d$tip))) Ntip(t2) + 1 else getMRCA(t2, tgt[[i]])
    lab <- paste0("FIX", i)
    t2 <- bind.tip(t2, lab, edge.length=0, where=nn)
    x[lab] <- as.character(fixnodes[[i]])
  }
  f <- fitMk(t2, x[t2$tip.label], fixedQ=Q, pi="fitzjohn")
  as.numeric(logLik(f))
}
res <- list(); sm_out <- list(); glm_out <- list()
for (v in c("any_S1", traits)) {
  y <- d[[v]]; names(y) <- d$tip
  cat("\n==== trait:", v, "====\n")
  tab <- table(d$cl, factor(y, levels=0:1)); print(tab)
  # dispersion within FSM: number of maximal all-carrier clades covering FSM carriers
  ft <- extract.clade(tt, nF); yf <- y[ft$tip.label]
  if (sum(yf) > 0) {
    nn <- (Ntip(ft)+1):(Ntip(ft)+ft$Nnode)
    pure <- sapply(nn, function(n) { tp <- extract.clade(ft,n)$tip.label; all(yf[tp]==1) })
    cover <- c(); for (n in nn[pure]) { tp <- extract.clade(ft,n)$tip.label; par <- ft$edge[ft$edge[,2]==n,1]; if (length(par)==0 || !all(yf[extract.clade(ft,par)$tip.label]==1)) cover <- c(cover, length(tp)) }
    singles <- sum(yf==1) - sum(cover)
    nblocks <- length(cover) + singles
    mono <- is.monophyletic(tt, names(y)[y==1 & d$cl=="FSM"])
    mrcaF <- if (sum(yf)>1) length(extract.clade(tt, getMRCA(tt, names(yf)[yf==1]))$tip.label) else 1
  } else { nblocks <- 0; mono <- NA; mrcaF <- NA; cover <- 0 }
  cat(sprintf("FSM carriers %d; monophyletic %s; MRCA of FSM carriers spans %s tips; maximal pure blocks %d (largest %d)\n",
      sum(yf), mono, mrcaF, nblocks, max(c(cover,1))))
  # parsimony
  p11 <- sankoff(y); p12 <- sankoff(y, gain=2, loss=1)
  mc <- mpr_counts(y); mc2 <- mpr_counts(y, wl=0.5)
  sc <- function(fx, g=1, l=1) sankoff(y, g, l, fx)
  F <- as.character(nF); FS <- as.character(nFS); R <- as.character(nR)
  scen <- list(stem_gain = setNames(list(1,0), c(F,FS)),
               later_gain = setNames(list(0), F),
               vertical_root = setNames(list(1), R),
               root_absent = setNames(list(0), R))
  pars <- sapply(scen, function(fx) c(eq=sc(fx), gain2=sc(fx,2,1)))
  cat("Sankoff cost (gain=loss=1):", p11, " (gain=2,loss=1):", p12, "\n"); print(pars)
  cat(sprintf("MPR (unordered) gains %d losses %d; MPR (loss cost 0.5) gains %d losses %d\n", mc$gains, mc$losses, mc2$gains, mc2$losses))
  # Mk
  fb <- fit_best(y); Q <- fb$best$index.matrix; Qm <- matrix(0,2,2); Qm[!is.na(Q)] <- fb$best$rates[Q[!is.na(Q)]]; diag(Qm) <- -rowSums(Qm)
  dimnames(Qm) <- list(fb$best$states, fb$best$states)
  cat("AIC ER/ARD:", round(fb$aic,2), " best", fb$bestname, " rates", signif(fb$best$rates,3), "\n")
  am <- asr_mk_model(tt, y+1, Nstates=2, transition_matrix=Qm[c("0","1"),c("0","1")], root_prior="stationary", include_ancestral_likelihoods=TRUE)
  pp <- function(n) am$ancestral_likelihoods[n-N, 2]
  cat(sprintf("P(present): ingroup root %.3f; FSM+sister node %.3f; FSM crown %.3f\n", pp(nR), pp(nFS), pp(nF)))
  L0 <- as.numeric(logLik(fb$best))
  cl_ <- sapply(scen, function(fx) tryCatch(cond_lik(y, Qm[c("0","1"),c("0","1")], fx), error=function(e) NA))
  cat("conditioned lnL:", round(cl_,2), " (unconditioned", round(L0,2), ")\n")
  # stochastic ancestral-state sampling from the joint posterior under the fitted Mk model (1,000 draws;
  # counts state changes between adjacent nodes, i.e. a lower bound on transitions per branch)
  cnt <- joint_sample_counts(y, Qm[c("0","1"),c("0","1")], 1000)
  q <- function(z) sprintf("%.1f [%g-%g]", mean(z), quantile(z,0.025), quantile(z,0.975))
  cat("joint posterior state samples (1000):", paste(colnames(cnt), apply(cnt,2,q), collapse="; "), "\n")
  # phyloglm, manuscript designs
  gfun <- function(sub) { dd <- d[d$cl %in% sub, ]; t2 <- drop.tip(tt, setdiff(tt$tip.label, dd$tip)); dd <- dd[t2$tip.label,]
    df <- data.frame(y=dd[[v]], FSM=as.integer(dd$cl=="FSM")); rownames(df) <- dd$tip
    f <- tryCatch(phyloglm(y ~ FSM, data=df, phy=t2, method="logistic_MPLE", btol=30), error=function(e) NULL)
    if (is.null(f)) return(c(n=nrow(dd), beta=NA, se=NA, p=NA, alpha=NA))
    s <- summary(f)$coefficients; c(n=nrow(dd), beta=s["FSM","Estimate"], se=s["FSM","StdErr"], p=s["FSM","p.value"], alpha=f$alpha) }
  g1 <- gfun(c("FSM","RGM","SGM")); g2 <- gfun(c("FSM","RGM"))
  cat(sprintf("phyloglm FSM vs RGM+SGM (n=%d): beta %.2f SE %.2f p %.3g alpha %.3g\n", g1["n"], g1["beta"], g1["se"], g1["p"], g1["alpha"]))
  cat(sprintf("phyloglm FSM vs RGM     (n=%d): beta %.2f SE %.2f p %.3g alpha %.3g\n", g2["n"], g2["beta"], g2["se"], g2["p"], g2["alpha"]))
  res[[v]] <- data.frame(trait=v, FSM_pos=tab["FSM","1"], FSM_n=sum(tab["FSM",]), RGM_pos=tab["RGM","1"], SGM_pos=tab["SGM","1"],
     FSM_carriers_monophyletic=mono, FSM_carrier_MRCA_tips=mrcaF, FSM_pure_blocks=nblocks, FSM_largest_block=max(c(cover,1)),
     sankoff_eq=p11, sankoff_gain2=p12, MPR_gains=mc$gains, MPR_losses=mc$losses,
     cost_stem_gain=pars["eq","stem_gain"], cost_later_gain=pars["eq","later_gain"], cost_vertical_root=pars["eq","vertical_root"], cost_root_absent=pars["eq","root_absent"],
     cost2_stem_gain=pars["gain2","stem_gain"], cost2_later_gain=pars["gain2","later_gain"], cost2_vertical_root=pars["gain2","vertical_root"], cost2_root_absent=pars["gain2","root_absent"],
     Mk_best=fb$bestname, AIC_ER=fb$aic["ER"], AIC_ARD=fb$aic["ARD"], q01=Qm["0","1"], q10=Qm["1","0"], lnL=L0,
     P_root=pp(nR), P_FSMsister=pp(nFS), P_FSMcrown=pp(nF),
     lnL_stem_gain=cl_["stem_gain"], lnL_later_gain=cl_["later_gain"], lnL_vertical_root=cl_["vertical_root"], lnL_root_absent=cl_["root_absent"],
     sm_gains=mean(cnt[,"gains"]), sm_losses=mean(cnt[,"losses"]), sm_gains_inFSM=mean(cnt[,"gains_inFSM"]), sm_losses_inFSM=mean(cnt[,"losses_inFSM"]),
     sm_gains_outFSM=mean(cnt[,"gains_outFSM"]), sm_frac_gain_on_FSM_stem=mean(cnt[,"gain_on_stem"]),
     glm_ingroup_n=g1["n"], glm_ingroup_beta=g1["beta"], glm_ingroup_se=g1["se"], glm_ingroup_p=g1["p"], glm_ingroup_alpha=g1["alpha"],
     glm_FSMRGM_n=g2["n"], glm_FSMRGM_beta=g2["beta"], glm_FSMRGM_se=g2["se"], glm_FSMRGM_p=g2["p"], glm_FSMRGM_alpha=g2["alpha"], row.names=NULL)
}
out <- do.call(rbind, res)
write.table(out, file.path(O, "species/lineage_species_tree_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
