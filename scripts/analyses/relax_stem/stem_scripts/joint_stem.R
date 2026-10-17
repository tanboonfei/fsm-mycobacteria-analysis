## Joint stem-edge loss probabilities for the 1,245 genus-conserved families (replaces the marginal product of asr178.R / asym3.R).
## Panels: 243 genomes (published) and 178 genomes >=95% complete. Models: per-family ER (castor, Table S6g) and the pooled-rate
## models of Table S6f (loss fixed at the pooled ER rate, gain = loss x {1, 1/2, 1/5, 1/10, 1/100, 2}).
suppressPackageStartupMessages({library(ape); library(castor); library(parallel)})
source("/path/to/project/analyses/relax_stem/stem/scripts/pruning.R")
B <- "/path/to/project/"; OUT <- paste0(B, "revision_checks/out/")
DST <- paste0(B, "analyses/relax_stem/stem/out/")
NC <- 2
mat <- read.delim(paste0(OUT, "m1_presence_matrix.tsv"), row.names = 1, check.names = FALSE)
meta <- read.delim(paste0(OUT, "m1_genome_meta.tsv"), stringsAsFactors = FALSE)
fates <- read.delim(paste0(OUT, "m1_family_fates.tsv"), stringsAsFactors = FALSE)
tr <- read.tree(paste0(B, "iqtree/myco_bac120.rooted.treefile"))
Qm <- function(g, l) matrix(c(-g, g, l, -l), 2, 2, byrow = TRUE)

panel <- function(minc) {             # identical construction to asr178.R
  keep <- meta$genome[meta$tip %in% tr$tip.label & meta$completeness >= minc]
  md <- meta[match(keep, meta$genome), ]; m <- mat[keep, , drop = FALSE]; rownames(m) <- md$tip
  tt <- drop.tip(tr, setdiff(tr$tip.label, md$tip)); m <- m[tt$tip.label, , drop = FALSE]
  md <- md[match(tt$tip.label, md$tip), ]
  crown <- getMRCA(tt, md$tip[md$clade == "FSM"]); parent <- tt$edge[tt$edge[, 2] == crown, 1]
  fsm_edges <- tt$edge[, 2][tt$edge[, 1] %in% c(crown, phytools::getDescendants(tt, crown))]   # child nodes of all edges inside the FSM crown subtree
  list(tt = tt, X = as.matrix(m), md = md, crown = crown, parent = parent, ntip = Ntip(tt), fsm_edges = fsm_edges,
       stem_len = tt$edge.length[tt$edge[, 2] == crown])
}
castor_marg <- function(P, Mfixed = NULL) {   # old method: castor rerooted marginals at crown and parent
  X <- P$X; nt <- P$ntip
  res <- mclapply(seq_len(ncol(X)), function(j) {
    x <- X[, j] + 1L
    if (length(unique(x)) == 1L) return(c(rep(as.numeric(unique(x) == 2L), 2), 0))
    f <- tryCatch(if (is.null(Mfixed)) asr_mk_model(P$tt, x, Nstates = 2, rate_model = "ER", include_ancestral_likelihoods = TRUE, reroot = TRUE)
                  else asr_mk_model(P$tt, x, Nstates = 2, transition_matrix = Mfixed, include_ancestral_likelihoods = TRUE), error = function(e) NULL)
    if (is.null(f) || !f$success) return(c(NA, NA, NA))
    c(f$ancestral_likelihoods[P$crown - nt, 2], f$ancestral_likelihoods[P$parent - nt, 2], f$transition_matrix[1, 2])
  }, mc.cores = NC)
  r <- do.call(rbind, res); colnames(r) <- c("p_crown", "p_parent", "q"); r
}
lost53 <- fates$family[fates$fate == "lost_cladewide"]
summ_rows <- list(); fam_rows <- list()
for (pn in list(list(minc = 0, tag = "243"), list(minc = 95, tag = "178"))) {
  P <- panel(pn$minc); X <- P$X; fam <- colnames(X); is53 <- fam %in% lost53
  cat(sprintf("\n=== panel %s: %d tips (FSM %d), crown %d, parent %d, stem length %.4f, %d families, FSM-subtree edges %d\n",
              pn$tag, P$ntip, sum(P$md$clade == "FSM"), P$crown, P$parent, P$stem_len, ncol(X), length(P$fsm_edges)))
  ## pooled ER rate (castor root_prior auto = max_likelihood, as in asym2.R)
  fe <- optimize(function(lr) -sum(mk_joint(P$tt, X, exp(lr), exp(lr))$loglik_ml), interval = c(log(1e-3), log(2)), tol = 1e-7)
  qE <- exp(fe$minimum); cat(sprintf("pooled ER rate %.4f (logL %.2f)\n", qE, -fe$objective))
  qE_use <- if (pn$tag == "243") 0.4610 else qE        # S6f used 0.4610 (243 panel); 178 panel uses its own pooled rate
  models <- list(
    list(name = "per-family ER (castor)", g = NA, l = NA),
    list(name = "pooled ER (gain = loss)", g = qE_use, l = qE_use),
    list(name = "gain = loss/2", g = qE_use / 2, l = qE_use),
    list(name = "gain = loss/5", g = qE_use / 5, l = qE_use),
    list(name = "gain = loss/10", g = qE_use / 10, l = qE_use),
    list(name = "gain = loss/100 (Dollo-like)", g = qE_use / 100, l = qE_use),
    list(name = "gain = loss x2", g = qE_use * 2, l = qE_use))
  for (mo in models) {
    t0 <- Sys.time()
    if (is.na(mo$g)) { cm <- castor_marg(P); g <- l <- cm[, "q"]; g[is.na(g)] <- 0; l <- g }
    else { cm <- castor_marg(P, Qm(mo$g, mo$l)); g <- rep(mo$g, ncol(X)); l <- rep(mo$l, ncol(X)) }
    old_prod <- cm[, "p_parent"] * (1 - cm[, "p_crown"])
    old95 <- cm[, "p_parent"] >= 0.95 & cm[, "p_crown"] <= 0.05; old90 <- cm[, "p_parent"] >= 0.90 & cm[, "p_crown"] <= 0.10
    out <- list()
    for (rp in c("stationary", "flat")) {
      r <- mk_joint(P$tt, X, g, l, root_prior = rp, edges = P$crown, event_edges = P$crown)
      J <- r$joint[[1]]; ev <- r$events
      rg <- mk_joint(P$tt, X, g, l, root_prior = rp, event_edges = P$fsm_edges)$events
      out[[rp]] <- list(J = J, ev = ev, gainFSM = rg[, "gain"], lossFSM = rg[, "loss"], mpar = r$marg[P$parent, ], mcr = r$marg[P$crown, ])
    }
    S <- out$stationary; jp <- S$J[, "p10"]
    na <- is.na(cm[, "p_crown"])
    summ_rows[[length(summ_rows) + 1]] <- data.frame(panel = pn$tag, model = mo$name,
      gain = if (is.na(mo$g)) "per family" else signif(mo$g, 4), loss = if (is.na(mo$l)) "per family" else signif(mo$l, 4),
      old_product_all = round(sum(old_prod, na.rm = TRUE), 1), old_conf95 = sum(old95, na.rm = TRUE), old_conf95_in53 = sum(old95 & is53, na.rm = TRUE),
      old_conf90 = sum(old90, na.rm = TRUE),
      joint_all = round(sum(jp, na.rm = TRUE), 1), joint_53 = round(sum(jp[is53], na.rm = TRUE), 1),
      joint_ge95 = sum(jp >= 0.95, na.rm = TRUE), joint_ge95_in53 = sum(jp >= 0.95 & is53, na.rm = TRUE),
      joint_ge90 = sum(jp >= 0.90, na.rm = TRUE), joint_ge90_in53 = sum(jp >= 0.90 & is53, na.rm = TRUE),
      Eloss_stem_all = round(sum(S$ev[, "loss"], na.rm = TRUE), 1), Eloss_stem_53 = round(sum(S$ev[is53, "loss"], na.rm = TRUE), 1),
      joint_all_flatroot = round(sum(out$flat$J[, "p10"], na.rm = TRUE), 1), joint_53_flatroot = round(sum(out$flat$J[is53, "p10"], na.rm = TRUE), 1),
      joint_ge95_flatroot = sum(out$flat$J[, "p10"] >= 0.95, na.rm = TRUE),
      castor_failures = sum(na), secs = round(as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    print(summ_rows[[length(summ_rows)]])
    fam_rows[[length(fam_rows) + 1]] <- data.frame(panel = pn$tag, model = mo$name, family = fam, fate = fates$fate[match(fam, fates$family)],
      prev_FSM_hc = fates$prev_FSM_hc[match(fam, fates$family)], gain = g, loss = l,
      castor_p_parent = round(cm[, "p_parent"], 5), castor_p_crown = round(cm[, "p_crown"], 5), old_product = round(old_prod, 5),
      old_conf95 = old95, post_p_parent = round(S$mpar, 5), post_p_crown = round(S$mcr, 5),
      joint_p10 = round(jp, 5), joint_p11 = round(S$J[, "p11"], 5), joint_p00 = round(S$J[, "p00"], 5), joint_p01 = round(S$J[, "p01"], 5),
      joint_p10_flatroot = round(out$flat$J[, "p10"], 5),
      E_loss_events_stem = round(S$ev[, "loss"], 5), E_gain_events_stem = round(S$ev[, "gain"], 5),
      E_gain_events_FSMcrown = round(S$gainFSM, 4), E_loss_events_FSMcrown = round(S$lossFSM, 4))
  }
  saveRDS(P, paste0(DST, "panel_", pn$tag, ".rds"))
}
summ <- do.call(rbind, summ_rows); famd <- do.call(rbind, fam_rows)
write.table(summ, paste0(DST, "joint_stem_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(famd, gzfile(paste0(DST, "joint_stem_per_family.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nwritten\n")
