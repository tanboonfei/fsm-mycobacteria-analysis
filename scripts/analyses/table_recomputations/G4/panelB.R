suppressPackageStartupMessages(library(ape))
P="/path/to/project/scripts/inputs/"
tr <- read.tree(paste0(P,"myco_bac120.rooted.treefile"))
d <- read.delim(paste0(P,"tsv/S1_genome_metadata.tsv"), stringsAsFactors=FALSE, na.strings=c("NA","","NaN","None"))
d$cl <- ifelse(grepl("^FSM",d$clade),"FSM", ifelse(grepl("^RGM",d$clade),"RGM", ifelse(d$clade=="SGM","SGM",NA)))
d <- d[!is.na(d$cl) & !duplicated(d$tip) & d$tip %in% tr$tip.label,]; rownames(d)<-d$tip
d$gs <- as.numeric(d$genome_size_Mb); d$comp <- as.numeric(d$completeness_used_pct)
ing <- drop.tip(tr, setdiff(tr$tip.label,d$tip))
run <- function(keep, label, thr=4){ x <- d[keep,]; t <- keep.tip(ing, x$tip); x <- x[t$tip.label,]
  sz <- x$gs/(x$comp/100); good <- x$cl!="FSM" & sz<thr; names(good)<-NULL
  pp <- prop.part(t); n <- length(t$tip.label)
  blocks <- list(); 
  # maximal blocks
  sets <- lapply(seq_along(pp), function(i) pp[[i]])
  allgood <- sapply(sets, function(s) all(good[s]))
  cand <- c(sets[allgood], as.list(which(good)))
  # maximal: not subset of another
  keepb <- sapply(seq_along(cand), function(i) !any(sapply(seq_along(cand), function(j) j!=i && length(cand[[j]])>length(cand[[i]]) && all(cand[[i]] %in% cand[[j]]))))
  mb <- cand[keepb]
  cat(label, "tips", n, "nonFSM small", sum(good), "blocks", length(mb), "largest", max(sapply(mb,length)), "\n") }
isMAG <- d$origin=="SPIRE MAG" | d$culture_status=="MAG / uncultured"
cat(table(d$origin),"\n")
ok <- !is.na(d$comp) & d$comp>=90
run(rownames(d)[ok], "B (all <90% pruned)")
okmag <- !is.na(d$comp) & !(d$comp<90 & d$culture_status=="MAG / uncultured")
run(rownames(d)[okmag], "B (only MAGs <90% pruned)")
cat("pruned MAGs", sum(d$comp<90 & d$culture_status=="MAG / uncultured",na.rm=TRUE), "pruned non-MAG", sum(d$comp<90 & d$culture_status!="MAG / uncultured",na.rm=TRUE), "NA comp", sum(is.na(d$comp)), "\n")
