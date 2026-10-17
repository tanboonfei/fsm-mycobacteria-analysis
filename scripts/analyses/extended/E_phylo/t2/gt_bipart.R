# frequency of competing bipartitions among the 299 gene trees (each target restricted to taxa present in that gene tree)
suppressMessages(library(ape))
gts <- read.tree("genes.treefile")
s1 <- read.delim("../t3/S1_genome_metadata.tsv", check.names=FALSE); cl <- setNames(s1$clade, s1$tip)
r1 <- read.tree("r1_lg.treefile"); tips <- r1$tip.label; og <- tips[grepl("Outgroup", cl[tips])]; r1 <- root(r1, og, resolve.root=TRUE)
fsm <- tips[grepl("^FSM", cl[tips])]; s42 <- readLines("sis_level1.txt"); pal <- readLines("sis_level2.txt")
m <- getMRCA(r1, s42); k <- r1$edge[r1$edge[,1]==m,2]; sub <- lapply(k, function(x) extract.clade(r1,x)$tip.label); sub <- sub[order(-lengths(sub))]
T <- list(`FSM`=fsm, `FSM + 42-genome clade (ML)`=c(fsm,s42), `42-genome clade`=s42, `FSM + rhodesiae subclade (H2)`=c(fsm,sub[[1]]), `FSM + sp020438785 subclade (H3)`=c(fsm,sub[[2]]),
          `FSM + palauense group (H4)`=c(fsm,pal), `42 + palauense (FSM outside; H5-like)`=c(s42,pal), `FSM + 42 + palauense`=c(fsm,s42,pal))
has <- function(tr, S) { S <- intersect(S, tr$tip.label); if (length(S)<2 || length(S)>=Ntip(tr)-1) return(NA)
  is.monophyletic(tr, S, reroot=TRUE) }
res <- sapply(T, function(S) sapply(gts, has, S=S))
out <- data.frame(bipartition=names(T), n_gene_trees_informative=colSums(!is.na(res)), n_present=colSums(res, na.rm=TRUE))
out$percent <- round(100*out$n_present/out$n_gene_trees_informative,1)
print(out, row.names=FALSE); write.table(out, "genetree_bipartition_freq.tsv", sep="\t", quote=FALSE, row.names=FALSE)
