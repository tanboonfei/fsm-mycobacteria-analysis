suppressPackageStartupMessages(library(ape))
R="/path/to/project/"
tr0=read.tree(paste0(R,"iqtree/myco_bac120.rooted.treefile"))
d=read.delim("S1.tsv",stringsAsFactors=FALSE); d=d[d$tip%in%tr0$tip.label,]
ing=d$tip[grepl("^FSM|^RGM",d$clade)|d$clade=="SGM"]; fsm=d$tip[grepl("^FSM",d$clade)]
tr=drop.tip(tr0,setdiff(tr0$tip.label,ing)); 
stemfrac=function(t){h=max(node.depth.edgelength(t)); n=getMRCA(t,fsm); e=which(t$edge[,2]==n); c(stem=t$edge.length[e],height=h,frac=t$edge.length[e]/h)}
print(stemfrac(tr))
tr$edge.length[tr$edge.length<=0]=1e-8
u=chronos(multi2di(tr),lambda=1,model="correlated")
print(stemfrac(u))
