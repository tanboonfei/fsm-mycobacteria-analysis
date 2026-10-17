# Builds the opsin tree object; support values are read by support_splits.R and opsin_splits.R (do not take support from this script after re-rooting).
suppressPackageStartupMessages(library(ape))
tr <- read.tree("../analysis_gaps/hgt/hgt_expanded_tree.contree")
s1 <- read.delim("../Publication/scripts/inputs/tsv/S1_genome_metadata.tsv", stringsAsFactors=FALSE)
cl <- setNames(ifelse(grepl("^FSM",s1$clade),"FSM",ifelse(grepl("^RGM",s1$clade),"RGM",ifelse(s1$clade=="SGM","SGM","other"))), sub("\\.fa$","",sub("^(GB_|RS_)","",s1$tip)))
lab <- tr$tip.label
gen <- sub("__.*$","", sub("\\|.*$","",lab)); gen <- sub("\\.fa$","",gen)
cls <- ifelse(gen %in% names(cl), cl[gen], NA)
acc <- sub("\\|.*$","",lab)
ref <- c(C0K2K6="actino",C0K2L2="actino",C0K2K8="actino",C0K2K7="actino",C0K2K9="actino",C0K2L3="actino",C0K2K5="actino",C0K2L0="actino",C0K2L1="actino",C0K2L4="actino",A0A249KS24="actino",A0A249L5J3="actino",
         Q6J4G7="PR",Q9F7P4="PR",Q9AFF7="PR",N0DLS5="PR",A4UHT3="PR",A4UHT6="PR",A0A075VQ34="PR",A4UHT4="PR",A4GIJ7="PR",A4GIH8="PR",A4UHT5="PR",F3L4J6="PR")
cls[is.na(cls) & acc %in% names(ref)] <- ref[acc[is.na(cls) & acc %in% names(ref)]]
cat("Tip classes:\n"); print(table(cls, useNA="ifany"))
cat("Unclassified tips:", head(lab[is.na(cls)], 15), "\n")
prs <- lab[grepl("Heliorhodopsin", lab)]
tr <- root(tr, outgroup=prs, resolve.root=TRUE)
fsm <- lab[cls %in% "FSM"]
m <- getMRCA(tr, fsm); inside <- extract.clade(tr, m)$tip.label
cat("\nFSM opsins:", length(fsm), "| MRCA clade size:", length(inside), "\n"); print(table(cls[match(inside, lab)], useNA="ifany"))
nl <- tr$node.label; cat("support at FSM-opsin MRCA:", nl[m-Ntip(tr)], "\n")
# largest pure-FSM clades
pure <- sapply((Ntip(tr)+1):(Ntip(tr)+tr$Nnode), function(n){ t <- extract.clade(tr,n)$tip.label; if (all(cls[match(t,lab)] %in% "FSM")) length(t) else 0 })
o <- order(pure, decreasing=TRUE)[1:5]; cat("largest FSM-only clades (size, support):", paste0(pure[o]," (",nl[o],")"), "\n")
big <- (Ntip(tr)+1:tr$Nnode)[o[1]]
par <- tr$edge[tr$edge[,2]==big,1]; sis <- setdiff(extract.clade(tr,par)$tip.label, extract.clade(tr,big)$tip.label)
cat("sister of largest FSM clade:", length(sis), "tips; classes:"); print(table(cls[match(sis,lab)], useNA="ifany")); cat("support of parent node:", nl[par-Ntip(tr)], "\n")
# actinorhodopsin clade relation
act <- lab[cls %in% "actino"]; ma <- getMRCA(tr, act); ia <- extract.clade(tr, ma)$tip.label
cat("\nMRCA of actinorhodopsin refs contains", length(ia), "tips:"); print(table(cls[match(ia,lab)], useNA="ifany"))
my <- lab[cls %in% c("FSM","RGM","SGM","other")]; mm <- getMRCA(tr, my); im <- extract.clade(tr, mm)$tip.label
cat("MRCA of all mycobacterial opsins contains", length(im), "tips:"); print(table(cls[match(im,lab)], useNA="ifany"))
write.tree(tr, "opsin_tree_rooted.nwk")
# group membership of the actinorhodopsin refs within FSM: are they a single clade nested among FSM?
pureA <- sapply((Ntip(tr)+1):(Ntip(tr)+tr$Nnode), function(n){ t <- extract.clade(tr,n)$tip.label; if (all(cls[match(t,lab)] %in% "actino")) length(t) else 0 })
cat("largest actino-only clade:", max(pureA), "of", length(act), "\n")
xpr <- lab[cls %in% "PR" | grepl("Xantho|Gloeobacter", lab)]
cat("Is (mycobacterial + actino) a clade excluding PR/xantho? ", !any(im %in% xpr), "\n")
