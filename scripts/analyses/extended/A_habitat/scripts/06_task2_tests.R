## Task 2: FSM vs RGM aquatic-biome and low-salinity associations, controlling for source study.
suppressPackageStartupMessages({library(lme4); library(logistf)})
set.seed(20261014)
O <- "/path/to/project/analyses/extended/A_habitat/"
d <- read.delim(paste0(O,"mag_master.tsv"), stringsAsFactors=FALSE, na.strings=c("","NA"))
d <- d[d$clade3 %in% c("FSM","RGM"),]
d$fsm <- as.integer(d$clade3=="FSM")
open_ocean <- c("TARA_Oceans_protists_metaG","Malaspina_marine")
res <- list()
fisher_row <- function(tab, lab, unit, n_units){
  f <- fisher.test(tab)
  data.frame(test=lab, unit=unit, FSM_pos=tab["1","1"], FSM_n=sum(tab["1",]), RGM_pos=tab["0","1"], RGM_n=sum(tab["0",]),
             OR=unname(f$estimate), CI_lo=f$conf.int[1], CI_hi=f$conf.int[2], p=f$p.value, n_units=n_units, note="")
}
mk <- function(x){ table(factor(x$fsm,levels=c(0,1)), factor(x$y,levels=c(0,1))) -> t; dimnames(t) <- list(c("0","1"),c("0","1")); t }
collapse_majority <- function(x, by){
  a <- aggregate(y ~ fsm + key, data=transform(x, key=by), FUN=mean)
  a$y <- as.integer(a$y >= 0.5); a }
run_all <- function(x, lab){
  out <- list()
  x$sp_study <- paste(x$study, x$species_key, sep="|")
  out[[1]] <- fisher_row(mk(x), lab, "genome", nrow(x))
  a <- collapse_majority(x, x$sp_study)
  out[[2]] <- fisher_row(mk(a), lab, "species x study (majority; ties->1)", nrow(a))
  # random one genome per species x study, 1000 draws
  ps <- c(); ors <- c()
  idx <- split(seq_len(nrow(x)), x$sp_study)
  for (i in 1:1000){ s <- x[sapply(idx, function(v) v[sample.int(length(v),1)]),]; f <- fisher.test(mk(s)); ps <- c(ps,f$p.value); ors <- c(ors,f$estimate) }
  r <- fisher_row(mk(x[sapply(idx, function(v) v[1]),]), lab, "species x study (random draw; 1000 reps)", length(idx))
  r$OR <- median(ors); r$CI_lo <- NA; r$CI_hi <- NA; r$p <- median(ps)
  r$note <- sprintf("median over draws; p range %.2g-%.2g; %% draws p<0.05 = %.0f", min(ps), max(ps), 100*mean(ps<0.05)); out[[3]] <- r
  b <- collapse_majority(x, x$study)
  out[[4]] <- fisher_row(mk(b), lab, "study x clade (majority)", nrow(b))
  # mixed model
  for (form in c("y ~ fsm + (1|study)", "y ~ fsm + (1|study) [nAGQ=20]", "y ~ fsm + (1|study) + (1|species_key)")){
    q <- if (grepl("nAGQ",form)) 20 else 1; ff <- sub(" \\[nAGQ=20\\]","",form)
    m <- tryCatch(glmer(as.formula(ff), data=x, family=binomial, nAGQ=q, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))), error=function(e) e)
    if (inherits(m,"error")) { out[[length(out)+1]] <- data.frame(test=lab, unit=paste("glmer", form), FSM_pos=NA,FSM_n=NA,RGM_pos=NA,RGM_n=NA,OR=NA,CI_lo=NA,CI_hi=NA,p=NA,n_units=nrow(x),note=conditionMessage(m)); next }
    m0 <- glmer(update(as.formula(ff), . ~ . - fsm), data=x, family=binomial, nAGQ=q, control=glmerControl(optimizer="bobyqa"))
    co <- summary(m)$coefficients["fsm",]
    lrt <- anova(m0, m)$`Pr(>Chisq)`[2]
    conv <- m@optinfo$conv$lme4$messages; sing <- isSingular(m)
    vc <- as.data.frame(VarCorr(m)); 
    degenerate <- max(vc$sdcor) > 10 || abs(co[1]) > 10
    note <- paste0(if (degenerate) "DEGENERATE (quasi-separation by study; estimates not interpretable); " else "", "Wald p=", signif(co[4],3), "; LRT p used; study SD=", round(vc$sdcor[vc$grp=="study"],2),
                   if (any(vc$grp=="species_key")) paste0("; species SD=", round(vc$sdcor[vc$grp=="species_key"],2)) else "",
                   if (sing) "; singular fit" else "", if (!is.null(conv)) paste0("; conv: ", paste(conv, collapse=" | ")) else "; converged")
    out[[length(out)+1]] <- data.frame(test=lab, unit=paste("glmer", form), FSM_pos=sum(x$y[x$fsm==1]), FSM_n=sum(x$fsm==1), RGM_pos=sum(x$y[x$fsm==0]), RGM_n=sum(x$fsm==0),
       OR=exp(co[1]), CI_lo=exp(co[1]-1.96*co[2]), CI_hi=exp(co[1]+1.96*co[2]), p=lrt, n_units=nrow(x), note=note)
  }
  # Firth logistic on study-collapsed units (handles separation)
  f <- logistf(y ~ fsm, data=b)
  out[[length(out)+1]] <- data.frame(test=lab, unit="Firth logistic, study x clade units", FSM_pos=sum(b$y[b$fsm==1]),FSM_n=sum(b$fsm==1),RGM_pos=sum(b$y[b$fsm==0]),RGM_n=sum(b$fsm==0),
     OR=exp(coef(f)[2]), CI_lo=exp(f$ci.lower[2]), CI_hi=exp(f$ci.upper[2]), p=f$prob[2], n_units=nrow(b), note="profile-likelihood CI")
  do.call(rbind,out)
}
## (1) aquatic biome (S1 'biome'), as in S12b
x <- d[!is.na(d$biome),]; x$y <- as.integer(x$biome=="aquatic")
res[[1]] <- run_all(x, "Aquatic biome (S1), all annotated")
xr <- x[!x$study %in% open_ocean,]
res[[2]] <- run_all(xr, "Aquatic biome (S1), excluding open-ocean studies (Tara, Malaspina)")
## (2) low salinity (S1 'salinity_habitat'), as in Fig. 4c
x <- d[!is.na(d$salinity_habitat),]; x$y <- as.integer(x$salinity_habitat=="low salinity")
res[[3]] <- run_all(x, "Low salinity (S1), all annotated")
xr <- x[!x$study %in% open_ocean,]
res[[4]] <- run_all(xr, "Low salinity (S1), excluding open-ocean studies (Tara)")
xr2 <- x[!x$study %in% c(open_ocean,"JGI_Salton_Sea_sediment"),]
res[[5]] <- run_all(xr2, "Low salinity (S1), excluding Tara and the Salton Sea hypersaline sediment study")
## (3) metadata-based aquatic class for all true MAGs (SPIRE + NCBI metagenome-derived), source-sample category from protist/A
cat3 <- read.delim("/path/to/project/scripts/analyses/protist/A/mag_source_metadata_summary.tsv", stringsAsFactors=FALSE)
d$srccat <- cat3$category[match(d$t, cat3$t)]
x <- d[d$source_set!="GTDB genome not flagged metagenome-derived" & !d$srccat %in% c("other/unknown"),]
x$y <- as.integer(x$srccat %in% c("water column","sediment"))
res[[6]] <- run_all(x, "Water column or sediment source (metadata category), all true MAGs")
xr <- x[!x$study %in% open_ocean,]
res[[7]] <- run_all(xr, "Water column or sediment source (metadata category), true MAGs excl. open ocean")
R <- do.call(rbind,res)
R$OR <- signif(R$OR,3); R$CI_lo <- signif(R$CI_lo,3); R$CI_hi <- signif(R$CI_hi,3); R$p <- signif(R$p,3)
write.table(R, paste0(O,"T2_clade_habitat_tests.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
options(width=250); print(R[,c("test","unit","FSM_pos","FSM_n","RGM_pos","RGM_n","OR","CI_lo","CI_hi","p","n_units")], row.names=FALSE)
cat("\nNotes:\n"); for (i in which(R$note!="")) cat(R$test[i]," | ",R$unit[i]," | ",R$note[i],"\n")
## studies shared by both clades among annotated sets
for (v in c("biome","salinity_habitat")){ z <- d[!is.na(d[[v]]),]; s <- tapply(z$clade3, z$study, function(q) length(unique(q))); cat(v,": studies with both clades:", paste(names(s)[s==2],collapse=", "),"\n") }
writeLines(capture.output(sessionInfo()), paste0(O,"logs/R_sessionInfo.txt"))
