# For control clades not monophyletic in the constrained 52-gene ML tree: build a full constraint tree by grafting the
# clade (bac120 internal topology) onto the primary 37-taxon topology at its bac120 position (sister group S restricted to panel taxa).
suppressMessages(library(ape)); library(jsonlite)
a <- commandArgs(TRUE); C <- a[1]; D <- "/path/to/project/analyses/extended/D_relax/"
x <- readRDS(paste0(D,"controls/clades.rds")); bt <- x$pt
cl <- fromJSON(paste0(D,"controls/control_clades.json"))[[C]]$tips
prim <- read.tree(paste0(D,"inputs/constraint.nwk")); p37 <- prim$tip.label
bt2 <- keep.tip(bt, c(intersect(p37, bt$tip.label), cl))
bt2 <- root(bt2, outgroup="spire_mag_00492514", resolve.root=FALSE)
stopifnot(is.monophyletic(bt2, cl))
m <- getMRCA(bt2, cl); par <- bt2$edge[bt2$edge[,2]==m,1]
sis <- setdiff(extract.clade(bt2, par)$tip.label, cl)
cat(C, "bac120 sister (panel taxa):", paste(sis, collapse=","), "\n")
sub <- keep.tip(bt2, cl); sub$edge.length <- NULL; sub$node.label <- NULL
pr <- root(prim, outgroup="spire_mag_00492514", resolve.root=FALSE)
if (length(sis) == 1) { where <- which(pr$tip.label == sis) } else { where <- getMRCA(pr, sis) }
if (where == Ntip(pr) + 1) stop("sister is whole tree")
pr$edge.length <- rep(1, nrow(pr$edge)); sub$edge.length <- rep(1, nrow(sub$edge)); sub$root.edge <- 1
g <- bind.tree(pr, sub, where = where, position = 0.5)
g$edge.length <- NULL; g$node.label <- NULL
stopifnot(is.monophyletic(g, cl))
write.tree(g, paste0(D, "controls/", C, ".constraint.nwk"))
