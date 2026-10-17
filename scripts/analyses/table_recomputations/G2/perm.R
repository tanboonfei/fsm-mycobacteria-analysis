suppressPackageStartupMessages({library(ape);library(phylolm)})
library(phangorn)
R="/path/to/project/"
tr0=read.tree(paste0(R,"iqtree/myco_bac120.rooted.treefile"))
d0=read.delim("S1.tsv",stringsAsFactors=FALSE); d0=d0[d0$tip %in% tr0$tip.label,]
d0$cl=ifelse(grepl("^FSM",d0$clade),"FSM",ifelse(grepl("^RGM",d0$clade),"RGM",ifelse(d0$clade=="SGM","SGM","other")))
d=d0[d0$cl%in%c("FSM","RGM","SGM"),]; tr=drop.tip(tr0,setdiff(tr0$tip.label,d$tip)); rownames(d)=d$tip; d=d[tr$tip.label,]
n=Ntip(tr); fsm=which(d$cl=="FSM"); 
# FSM node

pp=prop.part(tr); desc=c(as.list(1:n),pp)
sz=sapply((n+1):(n+tr$Nnode),function(i)length(desc[[i]]))
nodes=(n+1):(n+tr$Nnode)
cat("clade size distribution of internal nodes:\n"); print(table(cut(sz,c(0,5,10,20,50,67,100,134,200,268,400,823))))
fsmtips=fsm
fnode=nodes[sapply(nodes,function(i) setequal(desc[[i]],fsm))]; cat("FSM node",fnode,"\n")
vars=c(genome_size_Mb="size",coding_density_pct_uniform_prodigal="coding",GC_pct="GC",acidic_residues_pct="acidic")
for(v in names(vars)) d[[v]]=as.numeric(d[[v]])
res=list()
cand=nodes[sz>=67 & sz<=268]
cand=cand[!sapply(cand,function(i) any(desc[[i]]%in%fsm))|cand==fnode]
cat("candidate clades 67-268 tips not overlapping FSM (plus FSM):",length(cand),"\n")
for(i in cand) cat(i,sz[i-n],"RGMfrac",mean(d$cl[desc[[i]]]=="RGM"),"SGMfrac",mean(d$cl[desc[[i]]]=="SGM"),"size med",median(d$genome_size_Mb[desc[[i]]]),"\n")
# maximal (non-nested) candidate clades
maxc=cand[sapply(cand,function(i) !any(sapply(setdiff(cand,i),function(j) all(desc[[i]]%in%desc[[j]]) && j!=fnode && i!=fnode)))]
cat("maximal non-nested:",length(maxc),":",maxc,"\n")
stat=function(i,v){x=rep(0,n);x[desc[[i]]]=1;y=d[[v]];ok=!is.na(y);dd=data.frame(y=y[ok],x=x[ok],row.names=rownames(d)[ok]);tt=drop.tip(tr,setdiff(tr$tip.label,rownames(dd)))
 f=phylolm(y~x,data=dd,phy=tt,model="lambda");s=summary(f)$coefficients["x",];c(beta=s[1],t=s[3],p=s[4])}
out=NULL
for(v in names(vars)){ for(i in unique(c(fnode,cand))){ r=stat(i,v); out=rbind(out,data.frame(trait=vars[[v]],node=i,size=sz[i-n],beta=r[1],t=r[2],p=r[3]))}}
write.csv(out,"perm_clades.csv",row.names=FALSE)
for(tn in unique(out$trait)){o=out[out$trait==tn,];f=o[o$node==fnode,]
 for(set in list(all=unique(c(fnode,cand)),maximal=unique(c(fnode,maxc)))){oo=o[o$node%in%set,]
  cat(tn,"N=",nrow(oo),"FSM beta",round(f$beta,3),"t",round(f$t,2),"; p(|t|>=FSM)=",mean(abs(oo$t)>=abs(f$t)),"; p(beta<=FSM, one-sided)=",mean(oo$beta<=f$beta),"\n")}}
