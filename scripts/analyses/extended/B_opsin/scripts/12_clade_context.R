#!/usr/bin/env Rscript
# Walk from each mycobacterial clade towards the root; print what each successive ancestor adds, with bipartition support
suppressPackageStartupMessages(library(ape))
O <- "/path/to/project/analyses/extended/B_opsin"
u <- read.tree(file.path(O,"tree/opsin_ml.treefile")); r <- read.tree(file.path(O,"tree/opsin_ml.rooted_HeR.nwk"))
A <- read.delim(file.path(O,"tree/tip_annotation.tsv"), stringsAsFactors=FALSE); rownames(A) <- A$tip
Rz <- read.delim(file.path(O,"tree/tip_residues.tsv"), stringsAsFactors=FALSE); rownames(Rz) <- Rz$tip
M <- read.delim(file.path(O,"tree/mycobacterial_clades.tsv"), stringsAsFactors=FALSE)
lab <- u$tip.label; n <- Ntip(u)
key <- function(s) { s <- sort(s); c2 <- sort(setdiff(lab, s)); if (length(c2) < length(s) || (length(c2)==length(s) && c2[1] < s[1])) s <- c2; paste(s, collapse=",") }
sup <- list(); for (i in (n+1):(n+u$Nnode)) { s <- extract.clade(u, i)$tip.label; if (length(s) < n) sup[[key(s)]] <- u$node.label[i-n] }
getsup <- function(s) { k <- key(s); if (!is.null(sup[[k]])) sup[[k]] else "absent" }
nm <- function(t) ifelse(A[t,"set"]=="mycobacterial", paste0(A[t,"group"],":",t), paste0(t,"|",A[t,"label"],"|",A[t,"organism"],"|",A[t,"order"],"|",Rz[t,"motif"]))
N <- Ntip(r)
for (j in seq_len(nrow(M))) {
  mem <- strsplit(M$members[j], ",")[[1]]
  cur <- if (length(mem)>1) getMRCA(r, mem) else match(mem, r$tip.label)
  cat("\n######", M$clade_id[j], "n =", M$n[j], "\n"); prev <- mem
  for (step in 1:4) {
    p <- r$edge[r$edge[,2]==cur,1]; if (!length(p)) break
    s <- extract.clade(r,p)$tip.label; add <- setdiff(s, prev)
    cat(sprintf("-- ancestor %d: +%d tips, clade size %d, support %s\n", step, length(add), length(s), getsup(s)))
    mm <- add[A[add,"set"]=="mycobacterial"]; oo <- add[A[add,"set"]!="mycobacterial"]
    if (length(mm)) cat("   mycobacterial:", paste(names(table(A[mm,"group"])), table(A[mm,"group"]), collapse=" "), "\n")
    if (length(oo)) cat(paste0("   ", head(nm(oo), 40)), sep="\n")
    prev <- s; cur <- p
  }
}
