#!/usr/bin/env Rscript
# Lineage definition and support on the new opsin ML tree. Support is read from bipartitions of the
# unrooted IQ-TREE tree (never from re-rooted node labels). Rooting (for display/clade definition) on heliorhodopsins.
suppressPackageStartupMessages({library(ape)})
O <- "/path/to/project/analyses/extended/B_opsin"
args <- commandArgs(TRUE); tf <- if (length(args)) args[1] else file.path(O, "tree/opsin_ml.treefile")
u <- read.tree(tf)
A <- read.delim(file.path(O, "tree/tip_annotation.tsv"), stringsAsFactors=FALSE)
Rz <- read.delim(file.path(O, "tree/tip_residues.tsv"), stringsAsFactors=FALSE)
L <- read.delim(file.path(O, "screen/opsin_loci_residues.tsv"), stringsAsFactors=FALSE)
A <- merge(A, Rz[,c("tip","motif")], by="tip", all.x=TRUE); rownames(A) <- A$tip
lab <- u$tip.label; n <- Ntip(u)
# split -> support map from unrooted tree
key <- function(s) { s <- sort(s); c2 <- sort(setdiff(lab, s)); if (length(c2) < length(s) || (length(c2)==length(s) && c2[1] < s[1])) s <- c2; paste(s, collapse=",") }
sup <- list(); for (i in (n+1):(n+u$Nnode)) { s <- extract.clade(u, i)$tip.label; if (length(s) < n) sup[[key(s)]] <- u$node.label[i-n] }
getsup <- function(s) { k <- key(s); if (!is.null(sup[[k]])) sup[[k]] else "split absent" }
her <- A$tip[grepl("^HeR", A$label)]
r <- root(u, outgroup=her, resolve.root=TRUE); r$node.label <- NULL
cat("heliorhodopsins monophyletic in unrooted tree:", is.monophyletic(u, her), " support", getsup(her), "\n")
grp <- function(t) { a <- A[t,]; ifelse(a$set=="mycobacterial", paste0("Myco_", a$group), ifelse(a$set=="characterised", a$group, paste0(a$phylum,"/",a$order))) }
myco <- A$tip[A$set=="mycobacterial"]
# maximal mycobacteria-only clades
N <- Ntip(r); nodes <- (N+1):(N+r$Nnode)
pureM <- sapply(nodes, function(i) all(extract.clade(r,i)$tip.label %in% myco))
maxc <- c(); for (i in nodes[pureM]) { p <- r$edge[r$edge[,2]==i,1]; if (!(p %in% nodes[pureM])) maxc <- c(maxc, i) }
singles <- setdiff(myco, unlist(lapply(maxc, function(i) extract.clade(r,i)$tip.label)))
cat("maximal mycobacterial-only clades:", length(maxc), "; singleton mycobacterial tips:", length(singles), "\n")
desc_tbl <- function(t) { g <- table(grp(t)); paste(names(g), g, sep=":", collapse="; ") }
nearest_char <- function(i) { # walk up until a characterised reference is included
  cur <- i; repeat { p <- r$edge[r$edge[,2]==cur,1]; if (length(p)==0) return(NA); s <- extract.clade(r,p)$tip.label
    ch <- s[A[s,"set"]=="characterised"]; if (length(ch)) return(paste(A[ch,"label"], collapse=",")); cur <- p } }
rows <- list(); k <- 0
for (i in c(maxc, match(singles, r$tip.label))) {
  t <- if (i > N) extract.clade(r,i)$tip.label else r$tip.label[i]
  p <- r$edge[r$edge[,2]==i,1]; sis <- setdiff(extract.clade(r,p)$tip.label, t)
  k <- k+1
  rows[[k]] <- data.frame(clade_id=paste0("M", k), n=length(t), n_FSM=sum(A[t,"group"]=="FSM"), n_RGM=sum(A[t,"group"]=="RGM"), n_SGM=sum(A[t,"group"]=="SGM"),
    motifs=paste(names(table(A[t,"motif"])), table(A[t,"motif"]), sep=":", collapse=";"),
    support_clade=if (length(t)>1) getsup(t) else NA, sister_n=length(sis), sister_composition=desc_tbl(sis),
    support_clade_plus_sister=getsup(c(t,sis)), nearest_characterised=nearest_char(i), members=paste(t, collapse=","), stringsAsFactors=FALSE)
}
M <- do.call(rbind, rows); M <- M[order(-M$n),]; M$clade_id <- paste0("M", seq_len(nrow(M)))
write.table(M, file.path(O, "tree/mycobacterial_clades.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
print(M[, c("clade_id","n","n_FSM","n_RGM","n_SGM","motifs","support_clade","sister_n","support_clade_plus_sister","nearest_characterised")])
for (j in seq_len(min(nrow(M), 8))) cat(M$clade_id[j], "sister:", M$sister_composition[j], "\n")
write.tree(r, file.path(O, "tree/opsin_ml.rooted_HeR.nwk"))
