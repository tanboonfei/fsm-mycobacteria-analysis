# Recover the published lineage membership (65-sequence group vs 29 FSM-only clade) from the existing rooted opsin tree
suppressPackageStartupMessages(library(ape))
r <- read.tree("/path/to/project/scripts/analyses/opsin_tree_rooted.nwk")
lab <- r$tip.label
act <- lab[sub("\\|.*","",lab) %in% c("C0K2K6","C0K2L2","C0K2K8","C0K2K7","C0K2K9","C0K2L3","C0K2K5","C0K2L0","C0K2L1","C0K2L4","A0A249KS24","A0A249L5J3")]
pureA <- sapply((Ntip(r)+1):(Ntip(r)+r$Nnode), function(n){ t <- extract.clade(r,n)$tip.label; if (all(t %in% act)) length(t) else 0 })
na <- (Ntip(r)+1:r$Nnode)[which.max(pureA)]; a11 <- extract.clade(r,na)$tip.label; pa <- r$edge[r$edge[,2]==na,1]
s65 <- setdiff(extract.clade(r,pa)$tip.label, a11)
my <- lab[grepl("^spire_mag|^GB_|^RS_", lab)]
out <- data.frame(old_id=my, old_lineage=ifelse(my %in% s65, "L1_65group", "other"))
ia <- extract.clade(r, getMRCA(r, act))$tip.label
fsm_out <- my[!(my %in% ia)]
out$old_lineage[out$old_id %in% fsm_out] <- "L2_FSMonly"
write.table(out, "/path/to/project/analyses/extended/B_opsin/screen/old_tree_lineages.tsv", sep="\t", quote=FALSE, row.names=FALSE)
print(table(out$old_lineage))
