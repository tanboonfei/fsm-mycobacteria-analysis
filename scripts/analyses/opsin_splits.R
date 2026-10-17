suppressPackageStartupMessages(library(ape))
setwd("/path/to/project/analyses/")
L <- readLines("opsin_tree.R"); eval(parse(text=L[1:(grep("^tr <- root", L)-1)]))
u <- read.tree("../analysis_gaps/hgt/hgt_expanded_tree.treefile")
cat("first tip of contree:", u$tip.label[1], "\n")
look <- function(tr,tips){ n <- Ntip(tr)
  for (i in (n+1):(n+tr$Nnode)) { d <- extract.clade(tr,i)$tip.label
    if (setequal(d,tips) || setequal(setdiff(tr$tip.label,d),tips)) return(tr$node.label[i-n]) }
  "SPLIT ABSENT" }
r <- read.tree("opsin_tree_rooted.nwk")
act <- lab[cls %in% "actino"]; ma <- getMRCA(r, act); ia <- extract.clade(r, ma)$tip.label
pureA <- sapply((Ntip(r)+1):(Ntip(r)+r$Nnode), function(n){ t <- extract.clade(r,n)$tip.label; if (all(cls[match(t,lab)] %in% "actino")) length(t) else 0 })
na <- (Ntip(r)+1:r$Nnode)[which.max(pureA)]; a11 <- extract.clade(r,na)$tip.label; pa <- r$edge[r$edge[,2]==na,1]; s65 <- setdiff(extract.clade(r,pa)$tip.label, a11)
fout <- setdiff(lab[cls %in% "FSM"], ia)
cat("79-clade:", look(u, ia), "\n11 actino clade:", look(u, a11), "\n65 (60 FSM+5):", look(u, s65), "\n65+11:", look(u, c(s65,a11)), "\n29 FSM-only:", look(u, fout), "\n")
cat("heliorhodopsin tips:", sum(grepl("Heliorhodopsin", lab)), "\n")
# robustness to rooting: are these splits present irrespective of root? (they are bipartitions) -> check that each set is a clade under alternative roots
for (og in list(lab[grepl("Heliorhodopsin",lab)], lab[grepl("Channelrhodopsin",lab)], lab[grepl("Sensory|Halorhodopsin|BR_Halo",lab)])) {
  rr <- tryCatch(root(u, og, resolve.root=TRUE), error=function(e) NULL); if (is.null(rr)) {cat("root fail\n"); next}
  cat("root on", og[1], ": 65 clade", is.monophyletic(rr,s65), " 11", is.monophyletic(rr,a11), " 29", is.monophyletic(rr,fout), " 79", is.monophyletic(rr, ia), " 76", is.monophyletic(rr,c(s65,a11)), "\n") }
cat("non-FSM in 65:", s65[!cls[match(s65,lab)] %in% "FSM"], "\n")
