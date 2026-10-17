# read CF values for named bipartitions without re-rooting (labels stay on their nodes)
suppressMessages(library(ape))
a <- commandArgs(TRUE); pre <- a[1]
br <- read.tree(paste0(pre, ".cf.branch")); st <- read.table(paste0(pre, ".cf.stat"), header=TRUE, comment.char="#", check.names=FALSE)
s1 <- read.delim("../t3/S1_genome_metadata.tsv", check.names=FALSE); cl <- setNames(s1$clade, s1$tip); sp <- setNames(s1$species, s1$tip)
tips <- br$tip.label; fsm <- tips[grepl("^FSM", cl[tips])]
s42 <- intersect(readLines("sis_level1.txt"), tips); pal <- intersect(readLines("sis_level2.txt"), tips); sgm <- tips[cl[tips]=="SGM"]
r1 <- read.tree("r1_lg.treefile"); og <- r1$tip.label[grepl("Outgroup", cl[r1$tip.label])]; r1 <- root(r1, og, resolve.root=TRUE)
m <- getMRCA(r1, s42); k <- r1$edge[r1$edge[,1]==m,2]; sub <- lapply(k, function(x) extract.clade(r1,x)$tip.label); sub <- sub[order(-lengths(sub))]
targets <- list(FSM=fsm, `42-genome sister clade`=s42, `FSM + 42-genome clade`=c(fsm,s42), `FSM + 42 + palauense group`=c(fsm,s42,pal),
  `S42 subclade (rhodesiae, 30)`=sub[[1]], `S42 subclade (sp020438785/865, 12)`=sub[[2]], `FSM + rhodesiae subclade (H2)`=c(fsm,sub[[1]]), `FSM + sp020438785 subclade (H3)`=c(fsm,sub[[2]]), SGM=sgm)
desc <- function(n) if (n<=Ntip(br)) br$tip.label[n] else extract.clade(br,n)$tip.label
nn <- (Ntip(br)+1):(Ntip(br)+br$Nnode); dsets <- lapply(nn, desc)
out <- data.frame()
plyr_bind <- function(a,b){ if (!nrow(a)) return(b); for (n in setdiff(names(b),names(a))) a[[n]] <- NA; for (n in setdiff(names(a),names(b))) b[[n]] <- NA; rbind(a, b[,names(a)]) }
for (t in names(targets)) { T <- targets[[t]]; hit <- NA
  for (i in seq_along(nn)) { d <- dsets[[i]]; if (setequal(d,T) || setequal(setdiff(tips,d),T)) { hit <- br$node.label[i]; break } }
  row <- if (is.na(hit) || hit=="") NULL else st[st$ID==as.integer(hit),]
  out <- plyr_bind(out, data.frame(analysis=pre, bipartition=t, n_taxa=length(T), present_in_tree=!is.null(row) && nrow(row)>0,
    if (!is.null(row) && nrow(row)>0) row[,intersect(c("ID","gCF","gCF_N","gDF1","gDF2","gDFP","gN","sCF","sCF_N","sDF1","sDF2","sN","Label","Length"), names(row))] else data.frame(ID=NA)))
}
print(out, row.names=FALSE); write.table(out, paste0(pre, "_clades.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
