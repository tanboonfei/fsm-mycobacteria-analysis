# usage: Rscript label_control.R <iqtree treefile> <primary labelled nwk> <clade json> <clade name> <out nwk>
suppressMessages(library(ape)); library(jsonlite)
a <- commandArgs(TRUE)
rd <- function(f) { s <- readLines(f); s <- gsub("\\{([A-Za-z]+)\\}", "__\\1", s); read.tree(text=s) }
desc <- function(tr, n) if (n <= Ntip(tr)) tr$tip.label[n] else extract.clade(tr, n)$tip.label
prim <- rd(a[2])
# primary: map descendant-set key -> label (labels sit on child node / tip names)
plab <- list(); pnt <- Ntip(prim)
nm <- c(prim$tip.label, prim$node.label)
for (i in seq_len(nrow(prim$edge))) { ch <- prim$edge[i,2]
  lab <- sub(".*__", "", nm[ch]); if (!grepl("__", nm[ch])) lab <- "Unclassified"
  d <- sort(sub("__.*", "", desc(prim, ch))); plab[[paste(d, collapse=",")]] <- lab }
p37 <- sub("__.*", "", prim$tip.label)
tr <- read.tree(a[1]); tr$node.label <- NULL
tr <- root(tr, outgroup="spire_mag_00492514", resolve.root=FALSE)
cl <- fromJSON(a[3])[[a[4]]]$tips
stopifnot(is.monophyletic(tr, cl))
m <- getMRCA(tr, cl); nt <- Ntip(tr)
lab <- rep("", nrow(tr$edge))
for (i in seq_len(nrow(tr$edge))) { ch <- tr$edge[i,2]; d <- desc(tr, ch)
  if (all(d %in% cl)) { lab[i] <- if (ch <= nt || ch == m) "Test" else "Unclassified"; next }
  key <- paste(sort(intersect(d, p37)), collapse=",")
  l <- plab[[key]]; if (is.null(l)) stop(paste("no primary branch for", key))
  if (l == "Test") l <- "Unclassified"
  lab[i] <- l }
cat(a[4], "Test", sum(lab=="Test"), "Reference", sum(lab=="Reference"), "Unclassified", sum(lab=="Unclassified"),
    "test_len", round(sum(tr$edge.length[lab=="Test"]),5), "\n")
tl <- tr$tip.label; nl <- paste0("N", seq_len(tr$Nnode))
for (i in seq_len(nrow(tr$edge))) { ch <- tr$edge[i,2]
  if (ch <= nt) tl[ch] <- paste0(tl[ch], "{", lab[i], "}") else nl[ch-nt] <- paste0(nl[ch-nt], "{", lab[i], "}") }
tr$tip.label <- tl; tr$node.label <- nl
write.tree(tr, a[5])
