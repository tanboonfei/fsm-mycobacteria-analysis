#!/usr/bin/env Rscript
# Notung-like LCA reconciliation (duplication-loss only; D cost 1.5, L cost 1; gene tree rooted at the
# edge minimising cost) of each lineage's mycobacterial opsins against the 862-genome species tree.
suppressPackageStartupMessages(library(ape))
O <- "/path/to/project/analyses/extended/B_opsin"
P <- "/path/to/project/scripts"
sp <- read.tree(file.path(P,"inputs/myco_bac120.rooted.treefile")); sp$node.label <- NULL
u <- read.tree(file.path(O,"tree/opsin_ml.treefile"))
L <- read.delim(file.path(O,"screen/opsin_loci_lineages.tsv"), stringsAsFactors=FALSE); L <- L[L$in_tree %in% c("True","TRUE"),]
s1 <- read.delim(file.path(P,"inputs/tsv/S1_genome_metadata.tsv"), stringsAsFactors=FALSE); s1$key <- sub("\\.fa$","",s1$tip)
L$sptip <- s1$tip[match(L$genome, s1$key)]
fsmT <- s1$tip[grepl("^FSM", s1$clade)]
Ns <- Ntip(sp); par <- integer(Ns + sp$Nnode); par[sp$edge[,2]] <- sp$edge[,1]
dep <- node.depth.edgelength(sp); dep <- sapply(seq_along(par), function(i) { d <- 0; x <- i; while (par[x] != 0) { d <- d+1; x <- par[x] }; d })
anc <- function(x) { a <- x; while (par[x] != 0) { x <- par[x]; a <- c(a, x) }; a }
ancl <- lapply(seq_along(par), anc)
lca <- function(a, b) { A <- ancl[[a]]; for (x in ancl[[b]]) if (x %in% A) return(x) }
nF <- getMRCA(sp, fsmT)
recon <- function(g, map) {  # g rooted binary
  n <- Ntip(g); M <- integer(n + g$Nnode); M[1:n] <- match(map[g$tip.label], sp$tip.label)
  D <- 0; Ls <- 0; dupnodes <- c(); g <- reorder(g, "postorder")
  for (v in unique(g$edge[,1])) { ch <- g$edge[g$edge[,1]==v,2]; M[v] <- lca(M[ch[1]], M[ch[2]])
    isd <- M[v] %in% M[ch]; if (isd) { D <- D+1; dupnodes <- c(dupnodes, M[v]) }
    for (c in ch) Ls <- Ls + dep[M[c]] - dep[M[v]] - (if (isd) 0 else 1) }
  list(D=D, L=Ls, cost=1.5*D+Ls, root_map=M[n+1], dup_species_nodes=dupnodes)
}
best_root <- function(gu, map) {
  gu <- unroot(gu); gu$node.label <- NULL; best <- NULL
  for (e in seq_len(nrow(gu$edge))) { cn <- gu$edge[e,2]
    g <- if (cn <= Ntip(gu)) root(gu, outgroup=gu$tip.label[cn], resolve.root=TRUE) else root(gu, node=cn, resolve.root=TRUE)
    r <- recon(g, map); if (is.null(best) || r$cost < best$cost) best <- r }
  best
}
out <- list()
for (lin in c("A","B")) for (sub in c("all","FSM")) {
  x <- L[L$lineage==lin, ]; if (sub=="FSM") x <- x[x$clade=="FSM", ]
  g <- keep.tip(u, x$locus); map <- setNames(x$sptip, x$locus)
  b <- best_root(g, map)
  # reference: a species-tree-congruent gene tree for the same genomes (single copy each) -> D, L
  y <- x[!duplicated(x$sptip), ]; st <- keep.tip(sp, y$sptip); st$tip.label <- y$locus[match(st$tip.label, y$sptip)]
  ref <- recon(st, setNames(y$sptip, y$locus))
  rn <- b$root_map; rclade <- if (rn <= Ns) 1 else length(extract.clade(sp, rn)$tip.label)
  out[[paste(lin,sub)]] <- data.frame(lineage=lin, set=sub, genes=nrow(x), genomes=length(unique(x$sptip)),
     D=b$D, L=b$L, cost=b$cost, root_species_node_tips=rclade, root_is_FSM_crown=(rn==nF),
     root_within_FSM=(rn %in% c(nF, unlist(lapply(fsmT, function(t) ancl[[match(t, sp$tip.label)]])))) & all(extract.clade(sp, rn)$tip.label %in% fsmT) ,
     dups_within_FSM=sum(sapply(b$dup_species_nodes, function(d) d <= Ns && sp$tip.label[d] %in% fsmT || (d > Ns && all(extract.clade(sp,d)$tip.label %in% fsmT)))),
     congruent_ref_D=ref$D, congruent_ref_L=ref$L)
}
R <- do.call(rbind, out); print(R, row.names=FALSE)
write.table(R, file.path(O,"reconciliation/lca_dl_reconciliation.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
