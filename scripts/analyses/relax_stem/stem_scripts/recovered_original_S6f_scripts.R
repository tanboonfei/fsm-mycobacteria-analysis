## Original Table S6f scripts (not saved to disk at the time), recovered verbatim for provenance. Not re-run here.
## ===== recovered asym.R (run 2026-09-16T23:36:51.441Z; recovered copy of the script as run) =====
# Asymmetric (two-rate) gene-content model, with rates pooled across families so
# that they are identifiable, compared against the pooled equal-rates model.
suppressPackageStartupMessages({library(ape); library(castor); library(parallel)})
setwd("/path/to/project")
RE <- "analysis/reanalysis"
OUT <- "revision_checks/out"
set.seed(42)
m    <- read.delim(file.path(RE,"named_presence.tsv"), row.names=1, check.names=FALSE)
meta <- read.delim(file.path(RE,"named_meta.tsv"), stringsAsFactors=FALSE)
tr   <- read.tree("iqtree/myco_bac120.rooted.treefile")
keep <- meta$genome[meta$tip %in% tr$tip.label & meta$genome %in% rownames(m)]
md <- meta[match(keep, meta$genome),]
mm <- as.matrix(m[keep,,drop=FALSE]); rownames(mm) <- md$tip
tt <- drop.tip(tr, setdiff(tr$tip.label, md$tip))
mm <- mm[tt$tip.label,,drop=FALSE]; md <- md[match(tt$tip.label, md$tip),]
ntip <- length(tt$tip.label)
crown  <- getMRCA(tt, md$tip[md$clade=="FSM"])
parent <- tt$edge[tt$edge[,2]==crown,1]
cat(sprintf("tips=%d families=%d crown=%d parent=%d FSM=%d\n",
            ntip, ncol(mm), crown, parent, sum(md$clade=="FSM")))

states <- lapply(seq_len(ncol(mm)), function(j) mm[,j]+1L)
varying <- which(vapply(states, function(x) length(unique(x))>1L, logical(1)))
cat("families with variation:", length(varying), "\n")

Q <- function(q01,q10) matrix(c(-q01,q01,q10,-q10), 2, 2, byrow=TRUE)
NC <- 24
loglik <- function(par) {           # par = log(q01), log(q10)
  q01 <- exp(par[1]); q10 <- exp(par[2])
  if (!is.finite(q01) || !is.finite(q10) || q01>50 || q10>50) return(1e9)
  M <- Q(q01,q10)
  ll <- unlist(mclapply(states, function(x) {
         f <- tryCatch(asr_mk_model(tt, x, Nstates=2, transition_matrix=M,
                                    include_ancestral_likelihoods=FALSE),
                       error=function(e) NULL)
         if (is.null(f) || !f$success) return(NA_real_)
         f$loglikelihood }, mc.cores=NC))
  s <- sum(ll, na.rm=TRUE)
  if (!is.finite(s)) return(1e9)
  -s
}
cat("fitting pooled ARD (gain != loss) ...\n")
fitA <- optim(c(log(0.05), log(0.05)), loglik, method="Nelder-Mead",
              control=list(maxit=200, reltol=1e-8))
q01 <- exp(fitA$par[1]); q10 <- exp(fitA$par[2])
cat(sprintf("ARD: gain q01=%.5f  loss q10=%.5f  ratio loss/gain=%.2f  logL=%.2f\n",
            q01, q10, q10/q01, -fitA$value))
cat("fitting pooled ER (single rate) ...\n")
fitE <- optimize(function(lr) loglik(c(lr,lr)), interval=c(log(1e-4), log(5)), tol=1e-6)
qE <- exp(fitE$minimum)
cat(sprintf("ER : q=%.5f  logL=%.2f\n", qE, -fitE$objective))
lrt <- 2*((-fitA$value) - (-fitE$objective)); pv <- pchisq(lrt, df=1, lower.tail=FALSE)
cat(sprintf("LRT ARD vs ER = %.2f, df=1, p=%.3g\n", lrt, pv))

recon <- function(M, tag) {
  res <- mclapply(states, function(x) {
    if (length(unique(x))==1L) return(rep(as.numeric(unique(x)==2L),2))
    f <- tryCatch(asr_mk_model(tt, x, Nstates=2, transition_matrix=M,
                               include_ancestral_likelihoods=TRUE),
                  error=function(e) NULL)
    if (is.null(f) || !f$success) return(c(NA,NA))
    c(f$ancestral_likelihoods[crown-ntip,2], f$ancestral_likelihoods[parent-ntip,2])
  }, mc.cores=NC)
  P <- do.call(rbind,res); colnames(P) <- c("p_crown","p_preceding")
  cls <- function(p) ifelse(p>=0.95,"present",ifelse(p<=0.05,"absent","unresolved"))
  d <- data.frame(family=colnames(mm), p_crown=round(P[,1],4), p_preceding=round(P[,2],4),
                  crown=cls(P[,1]), preceding=cls(P[,2]), model=tag)
  cat("\n==", tag, "==\n"); print(table(preceding=d$preceding, crown=d$crown))
  cat(sprintf("confident stem-branch losses: %d\n",
              sum(d$preceding=="present" & d$crown=="absent", na.rm=TRUE)))
  cat(sprintf("reversals (absent -> present): %d\n",
              sum(d$preceding=="absent" & d$crown=="present", na.rm=TRUE)))
  cat(sprintf("unresolved at either node: %d\n", sum(d$crown=="unresolved"|d$preceding=="unresolved", na.rm=TRUE)))
  d
}
dA <- recon(Q(q01,q10), "pooled_ARD")
dE <- recon(Q(qE,qE),  "pooled_ER")
# Dollo-like: gain strongly penalised (1/100 of the fitted loss rate)
dD <- recon(Q(q10/100, q10), "gain_penalised_100x")
write.csv(rbind(dA,dE,dD), file.path(OUT,"asym_gene_content.csv"), row.names=FALSE)
saveRDS(list(q01=q01,q10=q10,qE=qE,lrt=lrt,p=pv), file.path(OUT,"asym_rates.rds"))
cat("\n-> revision_checks/out/asym_gene_content.csv\n")

## ===== recovered asym2.R (run 2026-09-16T23:38:33.151Z; recovered copy of the script as run) =====
suppressPackageStartupMessages({library(ape); library(castor); library(parallel)})
setwd("/path/to/project")
RE <- "analysis/reanalysis"; OUT <- "revision_checks/out"
set.seed(42)
m <- read.delim(file.path(RE,"named_presence.tsv"), row.names=1, check.names=FALSE)
meta <- read.delim(file.path(RE,"named_meta.tsv"), stringsAsFactors=FALSE)
tr <- read.tree("iqtree/myco_bac120.rooted.treefile")
keep <- meta$genome[meta$tip %in% tr$tip.label & meta$genome %in% rownames(m)]
md <- meta[match(keep, meta$genome),]
mm <- as.matrix(m[keep,,drop=FALSE]); rownames(mm) <- md$tip
tt <- drop.tip(tr, setdiff(tr$tip.label, md$tip)); mm <- mm[tt$tip.label,,drop=FALSE]
md <- md[match(tt$tip.label, md$tip),]; ntip <- length(tt$tip.label)
crown <- getMRCA(tt, md$tip[md$clade=="FSM"]); parent <- tt$edge[tt$edge[,2]==crown,1]
states <- lapply(seq_len(ncol(mm)), function(j) mm[,j]+1L)
Q <- function(g,l) matrix(c(-g,g,l,-l),2,2,byrow=TRUE)
NC <- 24
summ <- function(M, tag) {
  res <- mclapply(states, function(x){
    if (length(unique(x))==1L) return(rep(as.numeric(unique(x)==2L),2))
    f <- tryCatch(asr_mk_model(tt,x,Nstates=2,transition_matrix=M,
                  include_ancestral_likelihoods=TRUE), error=function(e) NULL)
    if (is.null(f)||!f$success) return(c(NA,NA))
    c(f$ancestral_likelihoods[crown-ntip,2], f$ancestral_likelihoods[parent-ntip,2])
  }, mc.cores=NC)
  P <- do.call(rbind,res)
  cl <- function(p) ifelse(p>=0.95,"present",ifelse(p<=0.05,"absent","unresolved"))
  cr <- cl(P[,1]); pr <- cl(P[,2])
  n_loss <- sum(pr=="present" & cr=="absent", na.rm=TRUE)
  n_rev  <- sum(pr=="absent"  & cr=="present", na.rm=TRUE)
  n_unr  <- sum(cr=="unresolved"|pr=="unresolved", na.rm=TRUE)
  cat(sprintf("%-26s stem_losses=%4d  reversals=%3d  unresolved=%4d\n", tag, n_loss, n_rev, n_unr))
  invisible(c(n_loss,n_rev,n_unr))
}
# pooled ER, refit with a sane upper bound
ll <- function(g,l){
  M<-Q(g,l)
  s<-sum(unlist(mclapply(states,function(x){
    f<-tryCatch(asr_mk_model(tt,x,Nstates=2,transition_matrix=M,
        include_ancestral_likelihoods=FALSE),error=function(e)NULL)
    if(is.null(f)||!f$success) NA_real_ else f$loglikelihood},mc.cores=NC)),na.rm=TRUE)
  if(!is.finite(s)) -1e12 else s }
fe <- optimize(function(lr) -ll(exp(lr),exp(lr)), interval=c(log(1e-3),log(2)), tol=1e-6)
qE <- exp(fe$minimum); cat(sprintf("pooled ER rate q = %.4f  logL = %.1f\n", qE, -fe$objective))
# asymmetric grid at fixed loss = qE, varying gain:loss ratio
cat("\nstem-branch losses as the gain rate is lowered relative to loss:\n")
res <- list()
for (r in c(1, 2, 5, 10, 100)) res[[as.character(r)]] <- summ(Q(qE/r, qE), sprintf("gain = loss/%d", r))
# and the reverse direction
for (r in c(2,5)) res[[paste0("hi",r)]] <- summ(Q(qE*r, qE), sprintf("gain = loss x%d", r))
saveRDS(list(qE=qE,res=res), file.path(OUT,"asym2.rds"))

## ===== recovered asym3.R (run 2026-09-16T23:40:50.890Z; recovered copy of the script as run) =====
suppressPackageStartupMessages({library(ape);library(castor);library(parallel)})
setwd("/path/to/project")
RE<-"analysis/reanalysis"; OUT<-"revision_checks/out"
m<-read.delim(file.path(RE,"named_presence.tsv"),row.names=1,check.names=FALSE)
meta<-read.delim(file.path(RE,"named_meta.tsv"),stringsAsFactors=FALSE)
tr<-read.tree("iqtree/myco_bac120.rooted.treefile")
keep<-meta$genome[meta$tip %in% tr$tip.label & meta$genome %in% rownames(m)]
md<-meta[match(keep,meta$genome),]; mm<-as.matrix(m[keep,,drop=FALSE]); rownames(mm)<-md$tip
tt<-drop.tip(tr,setdiff(tr$tip.label,md$tip)); mm<-mm[tt$tip.label,,drop=FALSE]
md<-md[match(tt$tip.label,md$tip),]; nt<-length(tt$tip.label)
cr<-getMRCA(tt,md$tip[md$clade=="FSM"]); pa<-tt$edge[tt$edge[,2]==cr,1]
states<-lapply(seq_len(ncol(mm)),function(j) mm[,j]+1L)
Q<-function(g,l) matrix(c(-g,g,l,-l),2,2,byrow=TRUE)
qE<-0.4610; NC<-24
row<-function(g,l,tag){
  P<-do.call(rbind,mclapply(states,function(x){
    if(length(unique(x))==1L){v<-as.numeric(unique(x)==2L);return(c(v,v))}
    f<-tryCatch(asr_mk_model(tt,x,Nstates=2,transition_matrix=Q(g,l),
        include_ancestral_likelihoods=TRUE),error=function(e)NULL)
    if(is.null(f)||!f$success) return(c(NA,NA))
    c(f$ancestral_likelihoods[cr-nt,2], f$ancestral_likelihoods[pa-nt,2])},mc.cores=NC))
  pc<-P[,1]; pp<-P[,2]
  exp_loss<-sum(pp*(1-pc),na.rm=TRUE)
  c95<-sum(pp>=0.95 & pc<=0.05,na.rm=TRUE); c90<-sum(pp>=0.90 & pc<=0.10,na.rm=TRUE)
  cat(sprintf("%-18s expected=%6.1f  n(0.95/0.05)=%4d  n(0.90/0.10)=%4d\n",tag,exp_loss,c95,c90))
  data.frame(model=tag,gain=g,loss=l,expected_stem_losses=round(exp_loss,1),
             n_conf_95=c95,n_conf_90=c90)
}
cat("Stem-branch losses of genus-conserved families under symmetric and asymmetric rates\n")
out<-rbind(
 row(qE,     qE, "ER (gain = loss)"),
 row(qE/2,   qE, "gain = loss/2"),
 row(qE/5,   qE, "gain = loss/5"),
 row(qE/10,  qE, "gain = loss/10"),
 row(qE/100, qE, "gain = loss/100 (Dollo-like)"),
 row(qE*2,   qE, "gain = loss x2"))
write.csv(out, file.path(OUT,"asym_stem_losses.csv"), row.names=FALSE)
cat("\n-> revision_checks/out/asym_stem_losses.csv\n")
