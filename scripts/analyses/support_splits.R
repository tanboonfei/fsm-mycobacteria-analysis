suppressMessages(library(ape))
a <- commandArgs(TRUE)
for (fn in a) {
 tr <- read.tree(fn); all <- tr$tip.label
 f <- intersect(readLines("fsm_tips.txt"), all); s <- intersect(readLines("sister_tips.txt"), all)
 desc <- prop.part(tr)  # per internal node, in node order
 lab <- tr$node.label
 getsup <- function(set){ for(i in seq_along(desc)){ d <- all[desc[[i]]]; if(setequal(d,set)||setequal(setdiff(all,d),set)) return(lab[i]) }; NA }
 cat(sprintf("%s | FSM: %s | FSM+sister: %s | sister alone: %s\n", basename(fn), getsup(f), getsup(union(f,s)), getsup(s)))
}
