#!/usr/bin/env Rscript
# Within each lineage on the rooted opsin tree: maximal FSM-only subclades, their support, and what interleaves them
suppressPackageStartupMessages(library(ape))
O <- "/path/to/project/analyses/extended/B_opsin"
u <- read.tree(file.path(O,"tree/opsin_ml.treefile")); r <- read.tree(file.path(O,"tree/opsin_ml.rooted_HeR.nwk"))
A <- read.delim(file.path(O,"tree/tip_annotation.tsv"), stringsAsFactors=FALSE); rownames(A) <- A$tip
M <- read.delim(file.path(O,"tree/mycobacterial_clades.tsv"), stringsAsFactors=FALSE)
lab <- u$tip.label; n <- Ntip(u)
key <- function(s) { s <- sort(s); c2 <- sort(setdiff(lab, s)); if (length(c2) < length(s) || (length(c2)==length(s) && c2[1] < s[1])) s <- c2; paste(s, collapse=",") }
sup <- list(); for (i in (n+1):(n+u$Nnode)) { s <- extract.clade(u, i)$tip.label; if (length(s) < n) sup[[key(s)]] <- u$node.label[i-n] }
getsup <- function(s) { k <- key(s); if (!is.null(sup[[k]])) sup[[k]] else "absent" }
mem <- function(id) strsplit(M$members[M$clade_id==id], ",")[[1]]
lin <- list(A = c(mem("M1"), mem("M3"), mem("M4")), B = mem("M2"))
N <- Ntip(r)
for (nm in names(lin)) {
  t <- lin[[nm]]; m <- getMRCA(r, t); cl <- extract.clade(r, m)
  cat(sprintf("\n== lineage %s: %d mycobacterial seqs; MRCA clade %d tips (support %s); non-mycobacterial inside: %d\n", nm, length(t), Ntip(cl), getsup(cl$tip.label), sum(!cl$tip.label %in% A$tip[A$set=="mycobacterial"])))
  isF <- A[cl$tip.label, "group"] == "FSM"; names(isF) <- cl$tip.label
  nn <- (Ntip(cl)+1):(Ntip(cl)+cl$Nnode)
  pure <- nn[sapply(nn, function(i) all(isF[extract.clade(cl,i)$tip.label]))]
  maxp <- pure[sapply(pure, function(i) { p <- cl$edge[cl$edge[,2]==i,1]; !(p %in% pure) })]
  sizes <- sapply(maxp, function(i) Ntip(extract.clade(cl,i)))
  fs <- sum(isF) - sum(sizes)
  cat(sprintf("FSM tips %d in %d maximal FSM-only clades (sizes %s) + %d singleton FSM tips\n", sum(isF), length(maxp), paste(sort(sizes, decreasing=TRUE), collapse=","), fs))
  for (i in maxp[order(-sizes)][1:min(6,length(maxp))]) { s <- extract.clade(cl,i)$tip.label; cat("  FSM clade n=", length(s), " support ", getsup(s), "\n") }
  cat("FSM opsins of lineage monophyletic:", is.monophyletic(r, t[A[t,"group"]=="FSM"]), "\n")
  write.tree(cl, file.path(O, paste0("tree/lineage", nm, "_subtree.nwk")))
}
