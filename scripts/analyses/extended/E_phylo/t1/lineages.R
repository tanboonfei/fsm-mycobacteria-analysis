suppressMessages(library(ape))
tr <- read.tree("Publication/scripts/analyses/reduced/r1_lg.treefile")
s1 <- read.delim("Publication/scripts/inputs/tsv/S1_genome_metadata.tsv", check.names=FALSE)
cl <- setNames(as.character(s1$clade), s1$tip); sp <- setNames(as.character(s1$species), s1$tip)
cla <- cl[tr$tip.label]; spp <- sp[tr$tip.label]
og <- tr$tip.label[grepl("Outgroup", cla)]
tr <- root(tr, og, resolve.root=TRUE)
cla <- cl[tr$tip.label]; spp <- sp[tr$tip.label]
fsm <- tr$tip.label[grepl("^FSM", cla)]
desc <- function(n) if (n<=Ntip(tr)) tr$tip.label[n] else extract.clade(tr,n)$tip.label
node <- getMRCA(tr, fsm); cat("FSM mono", length(desc(node))==length(fsm), "\n")
i <- 0
while (TRUE) {
  par <- tr$edge[tr$edge[,2]==node,1]
  sis <- setdiff(tr$edge[tr$edge[,1]==par,2], node)
  for (s in sis) { tt <- desc(s); i <- i+1
    cat(sprintf("level %d: n=%d clades=%s ; species e.g. %s\n", i, length(tt), paste(names(table(cla[tt])),table(cla[tt]),collapse=","), paste(head(unique(spp[tt]),6),collapse=", ")))
    writeLines(tt, sprintf("analyses/extended/E_phylo/t1/sis_level%d.txt", i)) }
  node <- par; if (node==Ntip(tr)+1) break
}
# also major sub-clades of the 42 sister clade
s42 <- readLines("Publication/scripts/analyses/sister_tips.txt"); s42 <- intersect(s42, tr$tip.label); cat("s42 in tree", length(s42), "\n")
m <- getMRCA(tr, s42); kids <- tr$edge[tr$edge[,1]==m,2]
for (k in kids) { tt <- desc(k); cat("S42 child n=", length(tt), paste(head(unique(spp[tt]),8),collapse=", "), "\n") }
