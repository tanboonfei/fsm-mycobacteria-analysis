## Task 3: opsin presence vs light/oxygen regime of the source sample, within FSM.
suppressPackageStartupMessages({library(logistf)})
set.seed(20261014)
O <- "/path/to/project/analyses/extended/A_habitat/"
x <- read.delim(paste0(O,"T3a_fsm_light_classification.tsv"), stringsAsFactors=FALSE)
x$op <- as.integer(x$rhodopsin=="present")
x$cls <- x$light_class_main
x$cls[grepl("^other",x$cls) | x$cls=="water: depth unresolved"] <- "other/unresolved"
lev <- c("photic water column","transition (10-20 m)","dark/anoxic water","sediment","wastewater/engineered","other/unresolved")
x$cls <- factor(x$cls, levels=lev)
collapse <- function(z){ z$u <- paste(z$study,z$species_key,z$cls,sep="|")
  a <- aggregate(op ~ u + cls + study + species_key, data=z, FUN=mean); a$op <- as.integer(a$op>=0.5); a }
summ <- function(z, unit){ t <- table(z$cls, factor(z$op,levels=0:1))
  data.frame(unit=unit, class=rownames(t), n=rowSums(t), opsin_pos=t[,2], pct=round(100*t[,2]/rowSums(t),1), n_studies=as.integer(tapply(z$study, z$cls, function(s) length(unique(s)))[rownames(t)])) }
x2 <- collapse(x)
S <- rbind(summ(x,"genome"), summ(x2,"species x study x class (majority)"))
contr <- list(
  "photic vs dark/anoxic water" = c("photic water column","dark/anoxic water"),
  "photic vs dark/anoxic+transition" = list("photic water column", c("dark/anoxic water","transition (10-20 m)")),
  "photic vs sediment" = c("photic water column","sediment"),
  "photic vs (dark/anoxic water + sediment)" = list("photic water column", c("dark/anoxic water","sediment")),
  "all natural aquatic (water+sediment) vs wastewater" = list(c("photic water column","transition (10-20 m)","dark/anoxic water","sediment"), "wastewater/engineered"))
T <- list()
for (dat in list(list(x,"genome"), list(x2,"species x study x class (majority)"))){
  z <- dat[[1]]
  zz <- droplevels(z[z$cls!="other/unresolved",])
  ft <- fisher.test(table(zz$cls, zz$op), simulate.p.value=FALSE, workspace=2e7)
  T[[length(T)+1]] <- data.frame(unit=dat[[2]], contrast="overall: 5 classes (excl. other/unresolved)", a_pos=NA,a_n=NA,b_pos=NA,b_n=NA, OR=NA, CI_lo=NA, CI_hi=NA, p=ft$p.value)
  zn <- droplevels(zz[zz$cls!="wastewater/engineered",])
  ft <- fisher.test(table(zn$cls, zn$op))
  T[[length(T)+1]] <- data.frame(unit=dat[[2]], contrast="overall: 4 natural aquatic classes (photic, transition, dark/anoxic, sediment)", a_pos=NA,a_n=NA,b_pos=NA,b_n=NA, OR=NA, CI_lo=NA, CI_hi=NA, p=ft$p.value)
  for (nm in names(contr)){ g <- contr[[nm]]; A <- z$cls %in% g[[1]]; B <- z$cls %in% g[[2]]
    tab <- matrix(c(sum(z$op[A]), sum(1-z$op[A]), sum(z$op[B]), sum(1-z$op[B])), 2, byrow=TRUE)
    f <- fisher.test(tab)
    T[[length(T)+1]] <- data.frame(unit=dat[[2]], contrast=nm, a_pos=tab[1,1], a_n=sum(tab[1,]), b_pos=tab[2,1], b_n=sum(tab[2,]), OR=unname(f$estimate), CI_lo=f$conf.int[1], CI_hi=f$conf.int[2], p=f$p.value) }
}
## sensitivity: photic restricted to documented depth/label (drop 'depth not reported'); transition merged into photic
xs <- x[x$light_class!="photic water column (depth not reported)",]
A <- xs$cls=="photic water column"; B <- xs$cls=="dark/anoxic water"; tab <- matrix(c(sum(xs$op[A]),sum(1-xs$op[A]),sum(xs$op[B]),sum(1-xs$op[B])),2,byrow=TRUE); f <- fisher.test(tab)
T[[length(T)+1]] <- data.frame(unit="genome", contrast="sensitivity: documented-photic only vs dark/anoxic", a_pos=tab[1,1],a_n=sum(tab[1,]),b_pos=tab[2,1],b_n=sum(tab[2,]),OR=unname(f$estimate),CI_lo=f$conf.int[1],CI_hi=f$conf.int[2],p=f$p.value)
A <- x$cls %in% c("photic water column","transition (10-20 m)"); B <- x$cls=="dark/anoxic water"; tab <- matrix(c(sum(x$op[A]),sum(1-x$op[A]),sum(x$op[B]),sum(1-x$op[B])),2,byrow=TRUE); f <- fisher.test(tab)
T[[length(T)+1]] <- data.frame(unit="genome", contrast="sensitivity: photic+transition vs dark/anoxic", a_pos=tab[1,1],a_n=sum(tab[1,]),b_pos=tab[2,1],b_n=sum(tab[2,]),OR=unname(f$estimate),CI_lo=f$conf.int[1],CI_hi=f$conf.int[2],p=f$p.value)
## Firth logistic, study-clustered units: opsin ~ dark (dark/anoxic+sediment vs photic), species x study units
zz <- x2[x2$cls %in% c("photic water column","dark/anoxic water","sediment"),]; zz$dark <- as.integer(zz$cls!="photic water column")
ff <- logistf(op ~ dark, data=zz)
T[[length(T)+1]] <- data.frame(unit="species x study x class (Firth)", contrast="dark/anoxic water + sediment vs photic", a_pos=sum(zz$op[zz$dark==1]),a_n=sum(zz$dark==1),b_pos=sum(zz$op[zz$dark==0]),b_n=sum(zz$dark==0),
   OR=exp(coef(ff)[2]),CI_lo=exp(ff$ci.lower[2]),CI_hi=exp(ff$ci.upper[2]),p=ff$prob[2])
T <- do.call(rbind,T); num <- c("OR","CI_lo","CI_hi","p"); T[num] <- lapply(T[num], signif, 3)
write.table(S, paste0(O,"T3b_opsin_by_light_class.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(T, paste0(O,"T3c_opsin_light_tests.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
options(width=220); print(S, row.names=FALSE); print(T, row.names=FALSE)
## species that occur in both photic and dark/sediment sources
sp <- tapply(as.character(x$cls), x$species_key, function(s) paste(sort(unique(s)),collapse=" + ")); print(sp[grepl("\\+",sp)])
