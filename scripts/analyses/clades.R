suppressMessages(library(ape))
tr <- read.tree("../Publication/scripts/inputs/myco_bac120.rooted.treefile")
s1 <- read.delim("../Publication/scripts/inputs/tsv/S1_genome_metadata.tsv", check.names=FALSE)
cat("tips", Ntip(tr), " S1 rows", nrow(s1), "\n")
cl <- setNames(as.character(s1$clade), s1$tip)
cat("tip match", sum(tr$tip.label %in% names(cl)), "\n")
fsm <- names(cl)[grepl("^FSM", cl)]
m <- getMRCA(tr, fsm); desc <- extract.clade(tr, m)$tip.label
cat("FSM n", length(fsm), " MRCA clade size", length(desc), " monophyletic", setequal(desc, fsm), "\n")
par <- tr$edge[tr$edge[,2]==m,1]
kids <- tr$edge[tr$edge[,1]==par,2]; sis <- setdiff(kids, m)
sisTips <- unlist(lapply(sis, function(k) if (k<=Ntip(tr)) tr$tip.label[k] else extract.clade(tr,k)$tip.label))
cat("sister n", length(sisTips), "\n")
d <- tr$edge[tr$edge[,1]==m,2]
dA <- if (d[1]<=Ntip(tr)) tr$tip.label[d[1]] else extract.clade(tr,d[1])$tip.label
dB <- setdiff(fsm, dA)
cat("FSM daughters", length(dA), length(dB), "\n")
nl <- tr$node.label; cat("support FSM node", nl[m-Ntip(tr)], " parent", nl[par-Ntip(tr)], "\n")
dn <- d[d>Ntip(tr)]; cat("daughter supports", nl[dn-Ntip(tr)], "\n")
writeLines(fsm, "fsm_tips.txt"); writeLines(sisTips, "sister_tips.txt"); writeLines(dA, "fsm_daughterA.txt"); writeLines(dB, "fsm_daughterB.txt")
writeLines(names(cl)[grepl("^RGM", cl)], "rgm_tips.txt"); writeLines(names(cl)[cl=="SGM"], "sgm_tips.txt")
writeLines(names(cl)[grepl("^Outgroup", cl)], "out_tips.txt"); writeLines(names(cl)[grepl("^Outside", cl)], "outside_tips.txt")
