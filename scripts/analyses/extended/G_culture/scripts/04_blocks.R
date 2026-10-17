# Small-genome block analysis (Table S13f, old S9f) under old and corrected culture labels.
# Block = maximal clade of the pruned ingroup tree whose tips are all small (< thr Mb) and of the
# same class (non-FSM or FSM). Reimplemented; the old-label run must reproduce the published S13f.
suppressMessages(library(ape))
G <- "/path/to/project/analyses/extended/G_culture/"
P <- "/path/to/project/scripts/inputs/"
s1 <- read.delim(paste0(P, "tsv/S1_genome_metadata.tsv"), check.names = FALSE)
cc <- read.delim(paste0(G, "corrected_culture_status.tsv"))
s1$new_status <- cc$new_status[match(s1$tip, cc$tip)]
s1$c3 <- sub(" .*", "", s1$clade)
newc <- read.delim(paste0(G, "new_completeness.tsv"))   # tip, completeness_used_new
s1$comp_new <- newc$completeness_used_new[match(s1$tip, newc$tip)]
tr <- read.tree(paste0(P, "myco_bac120.rooted.treefile"))
ing <- s1[s1$c3 %in% c("FSM", "RGM", "SGM"), ]
iso_old <- ing$culture_status %in% c("Cultured isolate", "NTM (isolate)", "Pathogen (isolate)") | ing$assembly_status == "closed"
iso_new <- grepl("^isolate", ing$new_status)

blocks <- function(tips, size, fsm, thr) {
  t <- keep.tip(tr, tips)
  sz <- size[match(t$tip.label, tips)]; f <- fsm[match(t$tip.label, tips)]
  res <- list()
  for (cls in c(FALSE, TRUE)) {
    small <- (sz < thr) & (f == cls)
    n <- Ntip(t); desc <- prop.part(t)   # tips under each internal node (node n+i)
    good <- sapply(desc, function(d) all(small[d]))
    # maximal: good node whose parent is not good
    par <- t$edge[match((n + 1):(n + t$Nnode), t$edge[, 2]), 1]
    gnode <- which(good & (is.na(par) | !good[par - n]))
    covered <- unique(unlist(desc[gnode]))
    singles <- setdiff(which(small), covered)
    D <- dist.nodes(t)
    sizes <- c(sapply(desc[gnode], length), rep(1, length(singles)))
    depth_max <- c(sapply(gnode, function(g) max(D[n + g, desc[[g]]])), rep(0, length(singles)))
    members <- c(lapply(desc[gnode], function(d) t$tip.label[d]), as.list(t$tip.label[singles]))
    res[[ifelse(cls, "FSM", "non")]] <- list(nsmall = sum(small), nblocks = length(sizes),
      largest = ifelse(length(sizes), max(sizes), 0), deep = sum(depth_max > 0.05), members = members,
      depth = depth_max, ntip = n, nFSM = sum(f))
  }
  res
}
panels <- function(lab, iso, comp) {
  out <- list(); mem <- list()
  A <- ing[iso, ]
  B <- ing[comp >= 90, ]; Bsize <- B$genome_size_Mb * 100 / comp[comp >= 90]
  for (thr in c(3.5, 4, 4.5)) {
    for (p in list(list("A: isolates and closed genomes, raw size", A$tip, A$genome_size_Mb, A$c3 == "FSM"),
                   list("B: all 824, completeness-corrected (MAGs <90% pruned)", B$tip, Bsize, B$c3 == "FSM"),
                   list("C: all 824, raw size (paper)", ing$tip, ing$genome_size_Mb, ing$c3 == "FSM"))) {
      r <- blocks(p[[2]], p[[3]], p[[4]], thr)
      out[[length(out) + 1]] <- data.frame(labels = lab, panel = p[[1]], threshold_Mb = thr, tips_in_tree = r$non$ntip,
        nonFSM_small_genomes = r$non$nsmall, nonFSM_blocks = r$non$nblocks, largest_nonFSM_block = r$non$largest,
        nonFSM_blocks_depth_gt_0.05 = r$non$deep, FSM_genomes_in_panel = r$non$nFSM,
        FSM_small_genomes = r$FSM$nsmall, largest_FSM_block = r$FSM$largest)
      if (thr == 4) mem[[p[[1]]]] <- r$non
    }
  }
  list(tab = do.call(rbind, out), mem = mem)
}
o <- panels("old (published S1)", iso_old, ing$completeness_used_pct)
n <- panels("corrected", iso_new, ing$comp_new)
tab <- rbind(o$tab, n$tab); tab <- tab[order(tab$panel, tab$threshold_Mb, tab$labels != "old (published S1)"), ]
write.table(tab, paste0(G, "S13f_blocks_old_vs_corrected.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(tab, row.names = FALSE)
# block membership at 4 Mb, both label sets
desc_tip <- function(x) { r <- ing[ing$tip == x, ]; sprintf("%s [%s; %s; %.2f Mb]", r$species, r$culture_status, r$new_status, r$genome_size_Mb) }
sink(paste0(G, "S13f_block_members_4Mb.txt"))
for (L in list(list("OLD", o), list("CORRECTED", n))) for (p in names(L[[2]]$mem)) {
  m <- L[[2]]$mem[[p]]; cat("\n==", L[[1]], p, ": ", m$nsmall, "small non-FSM genomes in", m$nblocks, "blocks\n")
  ord <- order(-sapply(m$members, length))
  for (i in ord) cat(length(m$members[[i]]), sprintf("depth %.4f", m$depth[i]), paste(sapply(m$members[[i]], desc_tip), collapse = "; "), "\n")
}
sink()
