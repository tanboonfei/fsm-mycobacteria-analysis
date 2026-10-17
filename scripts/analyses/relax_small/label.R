suppressMessages(library(ape)); library(jsonlite)
a <- commandArgs(TRUE); tr <- read.tree(a[1]); test <- fromJSON(a[2]); out <- a[3]
g <- fromJSON("/path/to/project/analysis_3way/groups.json")
ref <- c(g$RGM, g$SGM); lep <- c("leprae","lepromatosis","uberis")
tr$node.label <- NULL
tr <- root(tr, outgroup=lep, resolve.root=TRUE)
nt <- Ntip(tr); desc <- function(n) if (n<=nt) tr$tip.label[n] else extract.clade(tr,n)$tip.label
lab <- rep("", nrow(tr$edge)); m <- getMRCA(tr, test)
stopifnot(setequal(desc(m), test))
for (i in seq_len(nrow(tr$edge))) { ch <- tr$edge[i,2]; d <- desc(ch)
  if ((ch<=nt && d %in% test) || ch==m) lab[i] <- "{Test}"
  else if (all(d %in% ref)) lab[i] <- "{Reference}" }
# reference internal backbone: edges whose descendants are all ref already covered
cat("Test branches", sum(lab=="{Test}"), " Reference branches", sum(lab=="{Reference}"), "\n")
# write newick with labels appended to child names
nl <- rep("", tr$Nnode); tl <- tr$tip.label
for (i in seq_len(nrow(tr$edge))) { ch <- tr$edge[i,2]; if (lab[i]=="") next
  if (ch<=nt) tl[ch] <- paste0(tl[ch], lab[i]) else nl[ch-nt] <- lab[i] }
tr$tip.label <- tl; tr$node.label <- nl; tr$edge.length <- NULL
write.tree(tr, out)
