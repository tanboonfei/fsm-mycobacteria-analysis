## Validation of pruning.R on toy trees against brute-force enumeration, castor, and phytools stochastic mapping.
suppressPackageStartupMessages({library(ape); library(castor); library(phytools)})
source("/path/to/project/analyses/relax_stem/stem/scripts/pruning.R")
set.seed(1)
brute <- function(tree, x, g, l, pr1) {
  nt <- Ntip(tree); ni <- tree$Nnode; E <- tree$edge
  grid <- as.matrix(expand.grid(rep(list(0:1), ni)))
  st <- matrix(c(rep(x, each = nrow(grid))), nrow(grid)); st <- cbind(st, grid)   # columns = node numbers 1..nt+ni
  w <- ifelse(st[, nt + 1] == 1, pr1, 1 - pr1)
  for (k in seq_len(nrow(E))) { P <- Pmat(tree$edge.length[k], g, l)
    a <- st[, E[k, 1]]; b <- st[, E[k, 2]]
    w <- w * ifelse(a == 0 & b == 0, P$P00, ifelse(a == 0, P$P01, ifelse(b == 0, P$P10, P$P11))) }
  list(L = sum(w), st = st, w = w / sum(w))
}
maxdiff <- 0; ntests <- 0
for (rep in 1:30) {
  nt <- sample(5:8, 1); tree <- rtree(nt); if (rep %% 3 == 0) tree <- multi2di(tree)
  if (rep %% 5 == 0) tree <- di2multi(tree, tol = 0.3)        # include polytomies
  g <- runif(1, 0.05, 2); l <- runif(1, 0.05, 2); if (rep %% 2 == 0) g <- l   # ER and asymmetric
  x <- rbinom(nt, 1, 0.6)
  for (rp in c("stationary", "flat")) {
    pr1 <- if (rp == "flat") 0.5 else g / (g + l)
    b <- brute(tree, x, g, l, pr1)
    r <- mk_joint(tree, matrix(x, ncol = 1), g, l, root_prior = rp, edges = tree$edge[, 2])
    d1 <- abs(r$loglik - log(b$L))
    mb <- sapply(1:(nt + tree$Nnode), function(v) sum(b$w[b$st[, v] == 1]))
    d2 <- max(abs(r$marg[, 1] - mb))
    d3 <- 0
    for (k in seq_len(nrow(tree$edge))) { p <- tree$edge[k, 1]; c <- tree$edge[k, 2]
      J <- r$joint[[as.character(c)]]
      bj <- c(sum(b$w[b$st[, p] == 0 & b$st[, c] == 0]), sum(b$w[b$st[, p] == 0 & b$st[, c] == 1]),
              sum(b$w[b$st[, p] == 1 & b$st[, c] == 0]), sum(b$w[b$st[, p] == 1 & b$st[, c] == 1]))
      d3 <- max(d3, max(abs(J[1, ] - bj))) }
    maxdiff <- max(maxdiff, d1, d2, d3); ntests <- ntests + 1
  }
}
cat(sprintf("brute-force check: %d tree/pattern/prior combinations, max abs difference (loglik, node marginals, edge joints) = %.2e\n", ntests, maxdiff))

## castor equivalence (ER, where rerooted marginal likelihoods equal posteriors under the flat=stationary prior)
md <- 0; ml <- 0
for (rep in 1:20) { tree <- rtree(30); x <- rbinom(30, 1, 0.5); q <- runif(1, 0.1, 2)
  f <- asr_mk_model(tree, x + 1L, Nstates = 2, transition_matrix = matrix(c(-q, q, q, -q), 2, 2), include_ancestral_likelihoods = TRUE, reroot = TRUE)
  r <- mk_joint(tree, matrix(x, ncol = 1), q, q)
  md <- max(md, abs(f$ancestral_likelihoods[, 2] - r$marg[31:59, 1])); ml <- max(ml, abs(f$loglikelihood - r$loglik_ml)) }
cat(sprintf("castor ER check (20 random 30-tip trees): max abs difference in node marginals = %.2e; in loglik (castor root_prior auto = max_likelihood) = %.2e\n", md, ml))

## castor with an ASYMMETRIC matrix: rerooted 'likelihoods' are flat-prior conditional likelihoods, not posteriors
tree <- rtree(30); x <- rbinom(30, 1, 0.7); g <- 0.05; l <- 0.5
f <- asr_mk_model(tree, x + 1L, Nstates = 2, transition_matrix = matrix(c(-g, g, l, -l), 2, 2, byrow = TRUE), include_ancestral_likelihoods = TRUE, reroot = TRUE)
rs <- mk_joint(tree, matrix(x, ncol = 1), g, l, "stationary"); rf <- mk_joint(tree, matrix(x, ncol = 1), g, l, "flat")
cat(sprintf("asymmetric (gain=loss/10) 30-tip example: max |castor - stationary-prior posterior| = %.3f; castor vs flat-root posterior = %.3f\n",
  max(abs(f$ancestral_likelihoods[, 2] - rs$marg[31:59, 1])), max(abs(f$ancestral_likelihoods[, 2] - rf$marg[31:59, 1]))))

## expected loss events and joint posterior vs phytools stochastic mapping (ER and asymmetric), 12-tip tree
tree <- pbtree(n = 12, scale = 1); x <- setNames(c(1,1,1,0,0,1,0,1,1,1,0,1), tree$tip.label)
for (mod in list(c(g = 0.8, l = 0.8), c(g = 0.2, l = 1.0))) {
  g <- mod[["g"]]; l <- mod[["l"]]; pi1 <- g / (g + l)
  Q <- matrix(c(-g, g, l, -l), 2, 2, byrow = TRUE, dimnames = list(c("0","1"), c("0","1")))
  r <- mk_joint(tree, matrix(x, ncol = 1), g, l, "stationary", edges = tree$edge[, 2], event_edges = tree$edge[, 2])
  ## per-edge analytic loss events
  per_edge <- sapply(tree$edge[, 2], function(cn) mk_joint(tree, matrix(x, ncol = 1), g, l, "stationary", event_edges = cn)$events[1, "loss"])
  sm <- suppressMessages(make.simmap(tree, setNames(as.character(x), names(x)), Q = Q, pi = c(1 - pi1, pi1), nsim = 4000, message = FALSE))
  cnt <- t(sapply(sm, function(s) sapply(seq_len(nrow(s$edge)), function(k) { st <- names(s$maps[[k]]); sum(head(st, -1) == "1" & tail(st, -1) == "0") })))
  ## edge order of simmap trees equals input edge order
  sim_loss <- colMeans(cnt)
  jsim <- sapply(seq_len(nrow(tree$edge)), function(k) mean(sapply(sm, function(s) { st <- names(s$maps[[k]]); st[1] == "1" & tail(st, 1) == "0" })))
  jana <- sapply(tree$edge[, 2], function(cn) r$joint[[as.character(cn)]][1, "p10"])
  cat(sprintf("simmap check (gain %.1f, loss %.1f; 4000 maps): max |E[loss events] analytic - simmap| = %.3f (max analytic %.3f); max |P(1->0 endpoints) analytic - simmap| = %.3f\n",
              g, l, max(abs(per_edge - sim_loss)), max(per_edge), max(abs(jana - jsim))))
}
