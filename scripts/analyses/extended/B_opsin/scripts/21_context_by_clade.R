#!/usr/bin/env Rscript
# Do neighbourhood contexts (TetR/SDR site vs mannosidase/Rv0647c/AarF site vs Blh) sort by opsin subclade?
suppressPackageStartupMessages(library(ape))
O <- "/path/to/project/analyses/extended/B_opsin"
S <- read.delim(file.path(O,"synteny/opsin_neighbourhoods.tsv"), stringsAsFactors=FALSE); S[is.na(S)] <- ""
ctx <- apply(S[,c("g.2","g.1","g.1.1","g.2.1")], 1, function(v) { t <- paste(v, collapse=" ")
  if (grepl("carotene 15,15", t)) "Blh" else if (grepl("mannosidase|Rv0647c|AarF", t)) "MAN" else if (grepl("TetR|Short-chain|SDR", t)) "TET" else "other" })
names(ctx) <- S$locus
write.table(data.frame(locus=S$locus, context=ctx), file.path(O,"synteny/opsin_context.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
r <- read.tree(file.path(O,"tree/opsin_ml.rooted_HeR.nwk"))
for (lin in c("A","B")) {
  st <- read.tree(file.path(O, sprintf("tree/lineage%s_subtree.nwk", lin)))
  tips <- intersect(st$tip.label, names(ctx))
  cat("\nlineage", lin, "contexts:", paste(names(table(ctx[tips])), table(ctx[tips]), collapse=" "), "\n")
  for (k in c("TET","MAN","Blh")) { x <- tips[ctx[tips]==k]; if (length(x) < 2) next
    m <- getMRCA(st, x); cl <- extract.clade(st, m)$tip.label
    cat(sprintf("  %s: %d loci; monophyletic %s; MRCA clade %d tips containing contexts: %s\n", k, length(x), is.monophyletic(st, x), length(cl),
        paste(names(table(ctx[intersect(cl,names(ctx))])), table(ctx[intersect(cl,names(ctx))]), collapse=" "))) }
}
