suppressMessages(library(ape)); library(jsonlite)
B <- "/path/to/project/"
tr <- read.tree(paste0(B,"Publication/scripts/inputs/myco_bac120.rooted.treefile"))
s1 <- read.delim(paste0(B,"Publication/scripts/inputs/tsv/S1_genome_metadata.tsv"), check.names=FALSE, na.strings=c("NA","","NaN","None"))
s1$id <- sub("\\.fa$","",s1$tip)
genes <- sub("\\.ffn$","",list.files(paste0(B,"genes"), pattern="ffn$"))
g <- fromJSON(paste0(B,"analysis_3way/groups.json"))
panel <- c(g$RGM,g$SGM,g$FSM)
s1$panel <- s1$id %in% panel
s1$avail <- s1$id %in% genes & grepl("^(RGM|SGM)", s1$clade) & !s1$panel & s1$completeness_used_pct >= 70
cat("available non-panel RGM/SGM MAGs >=90% complete:", sum(s1$avail), "\n"); print(table(s1$clade[s1$avail]))
keep <- s1$tip[s1$avail | s1$panel]
cat("panel tips in tree:", sum(s1$tip[s1$panel] %in% tr$tip.label), "\n")
pt <- keep.tip(tr, intersect(keep, tr$tip.label))
lab <- setNames(s1$id, s1$tip); pt$tip.label <- lab[pt$tip.label]
isav <- setNames(s1$avail, s1$id)
nt <- Ntip(pt)
# candidate clades: internal nodes whose descendants are all available non-panel
cands <- list()
for (n in (nt+1):(nt+pt$Nnode)) { d <- extract.clade(pt, n)$tip.label
  if (all(isav[d]) && length(d) >= 3) cands[[length(cands)+1]] <- list(node=n, tips=d, n=length(d)) }
cat("candidate clades:", length(cands), "\n")
# maximal clades (not contained in a larger candidate)
sizes <- sapply(cands, `[[`, "n")
maxi <- Filter(function(c) !any(sapply(cands, function(o) o$n > c$n && all(c$tips %in% o$tips))), cands)
for (c in maxi) { ss <- s1[match(c$tips, s1$id),]
  cat(sprintf("node %d n=%d clade=%s species=%d size=%.2f-%.2f Mb compl=%.0f-%.0f\n", c$node, c$n, paste(unique(ss$clade),collapse="/"),
    length(unique(ss$species)), min(ss$genome_size_Mb), max(ss$genome_size_Mb), min(ss$completeness_used_pct), max(ss$completeness_used_pct))) }
saveRDS(list(pt=pt, maxi=maxi, cands=cands), "clades70.rds")
write.tree(pt, "pruned70.nwk")
