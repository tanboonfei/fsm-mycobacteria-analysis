## Two-state Mk pruning engine, vectorised across families (columns).
## State 0 = absent, 1 = present. gain = rate 0->1, loss = rate 1->0 (vectors of length F or scalars).
## Returns, for requested edges, the JOINT posterior of (parent state, child state) and expected numbers of
## gain / loss events on the edge, plus marginal posteriors at all nodes.
Pmat <- function(t, g, l) {             # returns list of 4 vectors P00,P01,P10,P11 (length F)
  s <- g + l; e <- exp(-s * t)
  pi1 <- ifelse(s > 0, g / s, 0.5); pi0 <- 1 - pi1
  e[s == 0] <- 1
  list(P00 = pi0 + pi1 * e, P01 = pi1 * (1 - e), P10 = pi0 * (1 - e), P11 = pi1 + pi0 * e)
}
## Expected number of k->l transitions on an edge of length t given start a and end b:
## E = q_kl * int_0^t P_ak(u) P_lb(t-u) du / P_ab(t); Gauss-Legendre quadrature (exact enough for smooth integrands)
gl <- local({ n <- 40; i <- 1:(n - 1); b <- i / sqrt(4 * i^2 - 1); J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- b; J[cbind(i + 1, i)] <- b; ev <- eigen(J, symmetric = TRUE)
  list(x = ev$values, w = 2 * ev$vectors[1, ]^2) })
expected_events <- function(t, g, l) {
  ## returns list: for each (a,b) in {00,01,10,11}: E[#loss], E[#gain]
  F <- length(g); out <- list()
  if (t <= 0) { z <- rep(0, F); for (ab in c("00","01","10","11")) out[[ab]] <- list(loss = z, gain = z); return(out) }
  u <- (gl$x + 1) * t / 2; w <- gl$w * t / 2
  Pt <- Pmat(t, g, l)
  get <- function(P, a, b) P[[paste0("P", a, b)]]
  for (a in 0:1) for (b in 0:1) {
    Iloss <- 0; Igain <- 0
    for (k in seq_along(u)) {
      Pu <- Pmat(u[k], g, l); Pr <- Pmat(t - u[k], g, l)
      Iloss <- Iloss + w[k] * get(Pu, a, 1) * get(Pr, 0, b)     # 1 -> 0 at time u
      Igain <- Igain + w[k] * get(Pu, a, 0) * get(Pr, 1, b)     # 0 -> 1 at time u
    }
    pab <- get(Pt, a, b)
    out[[paste0(a, b)]] <- list(loss = ifelse(pab > 0, l * Iloss / pab, 0), gain = ifelse(pab > 0, g * Igain / pab, 0))
  }
  out
}
## X: ntip x F 0/1 matrix with rows in tree$tip.label order. root_prior: "stationary" or "flat" (or numeric P(present)).
mk_joint <- function(tree, X, g, l, root_prior = "stationary", edges = NULL, event_edges = NULL) {
  F <- ncol(X); nt <- Ntip(tree); nn <- nt + tree$Nnode
  if (length(g) == 1) g <- rep(g, F); if (length(l) == 1) l <- rep(l, F)
  tr <- reorder(tree, "postorder"); E <- tr$edge; bl <- tr$edge.length
  D0 <- matrix(0, nn, F); D1 <- matrix(0, nn, F)
  D0[1:nt, ] <- (X == 0) * 1; D1[1:nt, ] <- (X == 1) * 1
  D0[(nt + 1):nn, ] <- 1; D1[(nt + 1):nn, ] <- 1
  logsc <- numeric(F)                 # log scale factors for the total likelihood
  M0 <- matrix(0, nrow(E), F); M1 <- matrix(0, nrow(E), F)   # messages child->parent per edge
  Pe <- vector("list", nrow(E))
  for (k in seq_len(nrow(E))) {      # postorder: children before parents
    p <- E[k, 1]; c <- E[k, 2]
    P <- Pmat(bl[k], g, l); Pe[[k]] <- P
    M0[k, ] <- P$P00 * D0[c, ] + P$P01 * D1[c, ]
    M1[k, ] <- P$P10 * D0[c, ] + P$P11 * D1[c, ]
    D0[p, ] <- D0[p, ] * M0[k, ]; D1[p, ] <- D1[p, ] * M1[k, ]
    # rescale parent when its last child has been processed (cheap: rescale every time)
    s <- D0[p, ] + D1[p, ]; s[s == 0] <- 1
    D0[p, ] <- D0[p, ] / s; D1[p, ] <- D1[p, ] / s; logsc <- logsc + log(s)
    M0[k, ] <- M0[k, ]; M1[k, ] <- M1[k, ]
  }
  ## Note: rescaling the parent after each child multiplies parent D by 1/s; product structure is preserved
  ## up to a per-family constant, which cancels in all posterior quantities.
  root <- nt + 1
  pr1 <- if (is.numeric(root_prior)) rep(root_prior, F) else if (root_prior == "flat") rep(0.5, F) else ifelse(g + l > 0, g / (g + l), 0.5)
  pr0 <- 1 - pr1
  lik_root <- pr0 * D0[root, ] + pr1 * D1[root, ]
  loglik <- log(lik_root) + logsc
  ## Upward pass: O_v(i) = P(data outside subtree(v), s_v = i), up to a per-family constant
  O0 <- matrix(0, nn, F); O1 <- matrix(0, nn, F); O0[root, ] <- pr0; O1[root, ] <- pr1
  ## sibling-excluded products S_p^{-c}(i) for each edge, computed directly from messages
  chl <- split(seq_len(nrow(E)), E[, 1])
  A0 <- matrix(0, nrow(E), F); A1 <- matrix(0, nrow(E), F)   # A_e(i) = O_p(i) * prod_{siblings} M_s(i)
  for (k in rev(seq_len(nrow(E)))) { # reverse postorder = preorder: parents before children
    p <- E[k, 1]; c <- E[k, 2]
    sib <- setdiff(chl[[as.character(p)]], k)
    a0 <- O0[p, ]; a1 <- O1[p, ]
    for (s in sib) { a0 <- a0 * M0[s, ]; a1 <- a1 * M1[s, ] }
    A0[k, ] <- a0; A1[k, ] <- a1
    P <- Pe[[k]]
    o0 <- a0 * P$P00 + a1 * P$P10; o1 <- a0 * P$P01 + a1 * P$P11
    sc <- o0 + o1; sc[sc == 0] <- 1
    O0[c, ] <- o0 / sc; O1[c, ] <- o1 / sc
  }
  marg1 <- O1 * D1; marg0 <- O0 * D0; marg <- marg1 / (marg0 + marg1)   # P(present | data) at each node
  joint_edge <- function(k) {
    c <- E[k, 2]; P <- Pe[[k]]
    j00 <- A0[k, ] * P$P00 * D0[c, ]; j01 <- A0[k, ] * P$P01 * D1[c, ]
    j10 <- A1[k, ] * P$P10 * D0[c, ]; j11 <- A1[k, ] * P$P11 * D1[c, ]
    tot <- j00 + j01 + j10 + j11
    cbind(p00 = j00 / tot, p01 = j01 / tot, p10 = j10 / tot, p11 = j11 / tot)
  }
  loglik_ml <- log(pmax(D0[root, ], D1[root, ])) + logsc   # castor root_prior='auto' (= 'max_likelihood')
  res <- list(loglik = loglik, loglik_ml = loglik_ml, marg = marg, edge = E, tree_po = tr)
  if (!is.null(edges)) {             # edges given as child node numbers (same numbering as input tree)
    res$joint <- lapply(edges, function(cn) joint_edge(which(E[, 2] == cn)))
    names(res$joint) <- edges
  }
  if (!is.null(event_edges)) {       # expected gain / loss events summed over the given child nodes' edges
    loss <- numeric(F); gain <- numeric(F)
    for (cn in event_edges) {
      k <- which(E[, 2] == cn); J <- joint_edge(k)
      ## expected_events requires per-family g,l; group families by identical rates for speed
      ee <- expected_events(bl[k], g, l)
      for (ab in c("00","01","10","11")) {
        w <- J[, paste0("p", ab)]
        loss <- loss + w * ee[[ab]]$loss; gain <- gain + w * ee[[ab]]$gain
      }
    }
    res$events <- cbind(loss = loss, gain = gain)
  }
  res
}
