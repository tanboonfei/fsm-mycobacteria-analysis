suppressMessages(library(ape))
s1 <- read.delim("Publication/scripts/inputs/tsv/S1_genome_metadata.tsv", check.names=FALSE)
cl <- setNames(as.character(s1$clade), s1$tip); sp <- setNames(as.character(s1$species), s1$tip)
d <- "analyses/extended/E_phylo/t1/"
for (f in c("Publication/scripts/analyses/reduced/r1_lg.treefile","Publication/scripts/analyses/reduced/r4_fast30.treefile","Publication/scripts/analyses/reduced/r4_fast30.contree","Publication/scripts/analyses/reduced/r3_pmsf_c20.treefile")) {
tr <- read.tree(f); og <- tr$tip.label[grepl("Outgroup", cl[tr$tip.label])]; tr <- root(tr, og, resolve.root=TRUE)
desc <- function(n) if (n<=Ntip(tr)) tr$tip.label[n] else extract.clade(tr,n)$tip.label
fsm <- tr$tip.label[grepl("^FSM", cl[tr$tip.label])]; node <- getMRCA(tr,fsm)
cat(f, "\n"); for (i in 1:3) { par <- tr$edge[tr$edge[,2]==node,1]; for (s in setdiff(tr$edge[tr$edge[,1]==par,2], node)) { tt <- desc(s); cat("  level",i,"n=",length(tt), paste(head(unique(sp[tt]),5),collapse=", "),"\n")}; node <- par }
}
tr <- read.tree("Publication/scripts/analyses/reduced/r1_lg.treefile"); og <- tr$tip.label[grepl("Outgroup", cl[tr$tip.label])]; tr <- root(tr, og, resolve.root=TRUE)
desc <- function(n) if (n<=Ntip(tr)) tr$tip.label[n] else extract.clade(tr,n)$tip.label
l3 <- readLines(paste0(d,"sis_level3.txt")); m <- getMRCA(tr,l3)
for (k in tr$edge[tr$edge[,1]==m,2]) { tt <- desc(k); cat("L3 child n=",length(tt), paste(head(unique(sp[tt]),10),collapse=", "),"\n") }
