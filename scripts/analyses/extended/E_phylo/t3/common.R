suppressPackageStartupMessages({library(ape)})
tr0 <- read.tree("myco_bac120.rooted.treefile")
d <- read.delim("S1_genome_metadata.tsv", stringsAsFactors=FALSE, check.names=FALSE)
d$cl <- ifelse(grepl("^FSM",d$clade),"FSM", ifelse(grepl("^RGM",d$clade),"RGM", ifelse(d$clade=="SGM","SGM",NA)))
rownames(d) <- d$tip
get_tree <- function(groups) {
  keep <- d$tip[d$cl %in% groups & !is.na(d$genome_size_Mb) & d$tip %in% tr0$tip.label]
  tt <- drop.tip(tr0, setdiff(tr0$tip.label, keep)); tt$node.label <- NULL
  tt <- reorder(tt, "postorder"); tt
}
fsm_tips <- d$tip[d$cl %in% "FSM"]
lep_tips <- d$tip[d$species %in% c("Mycobacterium leprae","Mycobacterium lepromatosis")]
s42 <- readLines("../../../../Publication/scripts/analyses/sister_tips.txt")
# describe the clade below an edge (edge index into tree$edge)
describe_edge <- function(tt, e, y) {
  n <- tt$edge[e,2]
  tips <- if (n <= Ntip(tt)) tt$tip.label[n] else extract.clade(tt, n)$tip.label
  cls <- table(d[tips,"cl"]); sp <- unique(d[tips,"species"])
  isFSM <- setequal(tips, intersect(fsm_tips, tt$tip.label))
  fsm_in <- sum(tips %in% fsm_tips)
  lab <- if (isFSM) "FSM stem" else if (all(tips %in% fsm_tips)) "within FSM" else
         if (length(intersect(lep_tips,tt$tip.label))>0 && all(intersect(lep_tips,tt$tip.label) %in% tips) && length(tips) <= 5) "leprosy clade (M. leprae, M. lepromatosis + relatives)" else
         if (setequal(tips, intersect(s42,tt$tip.label))) "42-genome sister clade" else if (setequal(tips, intersect(c(s42,fsm_tips),tt$tip.label))) "FSM + sister clade (parent of FSM stem)" else ""
  data.frame(edge=e, node=n, n_tips=length(tips), composition=paste(names(cls), cls, sep=":", collapse=","),
             FSM_tips_below=fsm_in, label=lab, mean_size_Mb_below=round(mean(d[tips,"genome_size_Mb"]),3),
             species_examples=paste(head(sp,4), collapse="; "), stringsAsFactors=FALSE)
}
# one genome per GTDB species (highest completeness, then fewest contigs); every "novel sp." genome kept
derep_tips <- function(groups) {
  x <- d[d$cl %in% groups & !is.na(d$genome_size_Mb),]
  x <- x[order(-as.numeric(x$completeness_used_pct), as.numeric(x$n_contigs)),]
  x$key <- ifelse(x$species=="novel sp.", x$tip, x$species)
  x$tip[!duplicated(x$key)]
}
