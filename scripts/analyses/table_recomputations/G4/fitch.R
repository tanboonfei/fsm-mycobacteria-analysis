suppressPackageStartupMessages({library(ape)})
P="/path/to/project/scripts/inputs/"
tr <- read.tree(paste0(P,"myco_bac120.rooted.treefile"))
d <- read.delim(paste0(P,"tsv/S1_genome_metadata.tsv"), stringsAsFactors=FALSE, na.strings=c("NA","","NaN","None"))
d$cl <- ifelse(grepl("^FSM",d$clade),"FSM", ifelse(grepl("^RGM",d$clade),"RGM", ifelse(d$clade=="SGM","SGM",NA)))
d <- d[!is.na(d$cl) & !duplicated(d$tip) & d$tip %in% tr$tip.label,]; rownames(d)<-d$tip
ing <- drop.tip(tr, setdiff(tr$tip.label, d$tip)); ing <- multi2di(ing)
st <- as.integer(as.numeric(d[ing$tip.label,"genome_size_Mb"])<4)
fitch <- function(tree, st){ n<-Ntip(tree); sets<-vector("list",n+tree$Nnode); for(i in 1:n) sets[[i]]<-st[i]
 ch <- split(tree$edge[,2], tree$edge[,1]); cost<-0
 for(nd in rev(unique(tree$edge[,1]))) {} 
 po <- postorder(tree); 
 for(e in seq(1,nrow(po),2)) { }
 NULL }
# simple recursive
cost <- 0
rec <- function(nd){ if(nd<=Ntip(ing)) return(st[nd]); k <- ing$edge[ing$edge[,1]==nd,2]; a<-rec(k[1]); b<-rec(k[2]); i<-intersect(a,b); if(length(i)) i else {cost<<-cost+1; union(a,b)} }
root <- Ntip(ing)+1; rec(root); cat("Fitch changes (<4 Mb vs >=4 Mb), all ingroup:", cost, "\n")
fsm <- d$tip[d$cl=="FSM"]; ft <- keep.tip(ing, fsm); stf <- as.integer(as.numeric(d[ft$tip.label,"genome_size_Mb"])<4)
cost2<-0
rec2 <- function(nd){ if(nd<=Ntip(ft)) return(stf[nd]); k <- ft$edge[ft$edge[,1]==nd,2]; a<-rec2(k[1]); b<-rec2(k[2]); i<-intersect(a,b); if(length(i)) i else {cost2<<-cost2+1; union(a,b)} }
rec2(Ntip(ft)+1); cat("Fitch changes within FSM subtree alone:", cost2, "\n")
