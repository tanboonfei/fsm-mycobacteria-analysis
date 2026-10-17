suppressMessages(library(ape))
tr <- ladderize(read.tree("../Publication/scripts/inputs/myco_bac120.rooted.treefile"))
ord <- tr$tip.label[tr$edge[tr$edge[,2] <= Ntip(tr), 2]]
pick <- function(file, k) { x <- readLines(file); o <- ord[ord %in% x]; if (length(o) <= k) return(o); o[unique(round(seq(1, length(o), length.out=k)))] }
B <- readLines("fsm_daughterB.txt")
fsm <- unique(c(B, pick("fsm_daughterA.txt", 35)))
sel <- unique(c(fsm, readLines("sister_tips.txt"), pick("rgm_tips.txt", 45), pick("sgm_tips.txt", 30), pick("outside_tips.txt", 4), pick("out_tips.txt", 4)))
cat("selected", length(sel), " FSM", length(fsm), "\n")
writeLines(sel, "reduced/taxa.txt"); writeLines(fsm, "reduced/fsm_sel.txt")
