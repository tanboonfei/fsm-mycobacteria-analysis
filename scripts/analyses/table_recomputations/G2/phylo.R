suppressPackageStartupMessages({library(ape);library(phylolm)})
R="/path/to/project/"
tr0=read.tree(paste0(R,"iqtree/myco_bac120.rooted.treefile"))
d0=read.delim("S1.tsv",stringsAsFactors=FALSE)
d0=d0[d0$tip %in% tr0$tip.label,]
d0$cl=ifelse(grepl("^FSM",d0$clade),"FSM",ifelse(grepl("^RGM",d0$clade),"RGM",ifelse(d0$clade=="SGM","SGM","other")))
sis=readLines(paste0(R,"analyses/sister_tips.txt"))
prep=function(keep){d=d0[keep,];tt=drop.tip(tr0,setdiff(tr0$tip.label,d$tip));rownames(d)=d$tip;d=d[tt$tip.label,];d$g=factor(ifelse(d$cl=="FSM","FSM","other"),levels=c("other","FSM"));list(d=d,t=tt)}
sets=list(n525=d0$cl%in%c("FSM","RGM"), n824=d0$cl%in%c("FSM","RGM","SGM"), n176=(d0$cl=="FSM")|(d0$tip %in% sis))
vars=c(genome_size_Mb="size",coding_density_pct_uniform_prodigal="coding",GC_pct="GC",acidic_residues_pct="acidic")
out=list()
for(nm in names(sets)){P=prep(sets[[nm]]);d=P$d;tr=P$t
 cat("==",nm,nrow(d),table(d$cl),"\n")
 for(v in names(vars)){y=suppressWarnings(as.numeric(d[[v]]));ok=!is.na(y);dd=d[ok,];dd$y=y[ok];tt=drop.tip(tr,setdiff(tr$tip.label,dd$tip))
  fit=phylolm(y~g,data=dd,phy=tt,model="lambda");s=summary(fit)$coefficients
  mw=wilcox.test(dd$y~dd$g)$p.value
  cat(sprintf("%s %s n=%d medFSM=%.4f medOth=%.4f MW=%.3g lambda=%.4f beta=%.4f se=%.4f t=%.3f p=%.4g\n",nm,vars[[v]],nrow(dd),median(dd$y[dd$g=="FSM"]),median(dd$y[dd$g=="other"]),mw,fit$optpar,s["gFSM",1],s["gFSM",2],s["gFSM",3],s["gFSM",4]))}
 for(v in c("rhodopsin","urease","ectoine_synthase")){yb=as.integer(d[[v]]=="present");dd=d;dd$y=yb
  tab=table(dd$g,dd$y);fp=fisher.test(tab)$p.value
  f=phyloglm(y~g,data=dd,phy=tr,method="logistic_MPLE",btol=30);cf=summary(f)$coefficients
  cat(sprintf("%s %s FSM %d/%d (%.1f%%) other %d/%d (%.1f%%) Fisher=%.3g beta=%.4f se=%.4f p=%.4g alpha=%.4f\n",nm,v,sum(yb[d$g=="FSM"]),sum(d$g=="FSM"),100*mean(yb[d$g=="FSM"]),sum(yb[d$g!="FSM"]),sum(d$g!="FSM"),100*mean(yb[d$g!="FSM"]),fp,cf["gFSM",1],cf["gFSM",2],cf["gFSM",4],f$alpha))}
}
