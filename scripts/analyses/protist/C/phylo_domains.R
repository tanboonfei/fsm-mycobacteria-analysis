suppressMessages({library(ape);library(phylolm)})
t<-read.tree("../../../Publication/scripts/inputs/myco_bac120.rooted.treefile")
d<-read.delim("domain_counts.tsv"); d<-d[d$group %in% c("FSM","RGM_sister","RGM_other") & d$comp>=90,]
rownames(d)<-d$tip; d<-d[d$tip %in% t$tip.label,]
t<-keep.tip(t,d$tip); d<-d[t$tip.label,]
d$FSM<-as.integer(d$group=="FSM"); d$lsize<-log(d$size)
res<-NULL
for(c in c("ELP_any","SET","PE","PPE","MCE","WXG100","Patatin","PLC")){
  d$y<-log1p(d[[c]])
  f<-phylolm(y~FSM+lsize+comp,data=d,phy=t,model="lambda")
  s<-summary(f)$coefficients
  f0<-lm(y~FSM+lsize+comp,data=d); s0<-summary(f0)$coefficients
  res<-rbind(res,data.frame(domain=c,n=nrow(d),FSM_coef_phylo=s["FSM",1],p_phylo=s["FSM",4],lambda=f$optpar,FSM_coef_ols=s0["FSM",1],p_ols=s0["FSM",4],
    FSM_median=median(d[[c]][d$FSM==1]),RGM_median=median(d[[c]][d$FSM==0])))
}
print(res,digits=3); write.table(res,"domain_phylolm.tsv",sep="\t",quote=F,row.names=F)
