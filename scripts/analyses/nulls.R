suppressPackageStartupMessages({library(ape);library(phytools);library(geiger);library(phylolm)})
tr <- read.tree("../Publication/scripts/inputs/myco_bac120.rooted.treefile")
d  <- read.delim("../Publication/scripts/inputs/tsv/S1_genome_metadata.tsv", stringsAsFactors=FALSE, na.strings=c("NA","","NaN","None"))
d$cl <- ifelse(grepl("^FSM",d$clade),"FSM", ifelse(grepl("^RGM",d$clade),"RGM", ifelse(d$clade=="SGM","SGM",NA)))
d$gs <- as.numeric(d$genome_size_Mb)
sink("nulls_out.txt", split=TRUE)
## 1. Brownian null on log genome size
set.seed(42)
bm <- function(groups, tag, transform) {
  dd <- d[d$cl %in% groups & !is.na(d$gs) & d$tip %in% tr$tip.label,]
  tt <- drop.tip(tr, setdiff(tr$tip.label, dd$tip)); rownames(dd) <- dd$tip; dd <- dd[tt$tip.label,]
  y <- setNames(transform(dd$gs), dd$tip); fsm <- dd$tip[dd$cl=="FSM"]; oth <- setdiff(names(y), fsm)
  obs <- mean(y[fsm]) - mean(y[oth])
  fit <- fitContinuous(tt, y, model="BM"); sims <- fastBM(tt, sig2=fit$opt$sigsq, nsim=2000, a=mean(y))
  nul <- apply(sims, 2, function(z) mean(z[fsm]) - mean(z[oth]))
  p <- (sum(abs(nul) >= abs(obs))+1)/2001
  rng <- range(sims); cat(sprintf("BM null %-20s %-6s n=%d obs %+.3f null sd %.3f p = %.3f | simulated tip range %.2f to %.2f\n", tag, deparse(substitute(transform)), length(y), obs, sd(nul), p, rng[1], rng[2]))
}
bm(c("FSM","RGM"), "FSM vs RGM", identity); bm(c("FSM","RGM"), "FSM vs RGM", log)
bm(c("FSM","RGM","SGM"), "FSM vs ingroup", identity); bm(c("FSM","RGM","SGM"), "FSM vs ingroup", log)
## 2. Opsin null conditioned on prevalence
set.seed(20260901)
d$opsin <- ifelse(is.na(d$rhodopsin), NA, as.integer(d$rhodopsin=="present"))
or_of <- function(cl, y){ tab <- table(factor(cl,levels=c("RGM","FSM")), factor(y,levels=c(0,1)))+0.5; (tab[2,2]*tab[1,1])/(tab[2,1]*tab[1,2]) }
dd <- d[d$cl %in% c("FSM","RGM") & !is.na(d$opsin) & d$tip %in% tr$tip.label,]
tt <- drop.tip(tr, setdiff(tr$tip.label, dd$tip)); rownames(dd) <- dd$tip; dd <- dd[tt$tip.label,]
y <- setNames(dd$opsin, dd$tip); npos <- sum(y); obs <- or_of(dd$cl, y)
rate <- fitMk(tt, setNames(as.character(y), names(y)), model="ER")$rates[1]
cat(sprintf("\nOpsin: n=%d positives=%d observed OR=%.2f fitted ER rate=%.4f\n", length(y), npos, obs, rate))
# original design: fraction monomorphic / median OR
orig <- replicate(2000, { s <- as.integer(as.character(rTraitDisc(tt, model="ER", k=2, rate=rate, states=c("0","1"), root.value=sample(1:2,1))[tt$tip.label])); c(sum(s), or_of(dd$cl,s)) })
cat(sprintf("Original design: replicates with <=5 or >=%d positives (near-monomorphic) = %.1f%%; median OR %.2f\n", length(y)-5, 100*mean(orig[1,]<=5 | orig[1,]>=length(y)-5), median(orig[2,])))
# conditioned: rejection sampling across rate multipliers, accept positives within +-10% of observed
acc <- c(); tried <- 0; lo <- floor(npos*0.9); hi <- ceiling(npos*1.1)
mults <- c(1,2,5,10,20,50)
for (m in mults) for (i in 1:20000) { tried <- tried+1
  s <- as.integer(as.character(rTraitDisc(tt, model="ER", k=2, rate=rate*m, states=c("0","1"), root.value=sample(1:2,1))[tt$tip.label]))
  if (sum(s)>=lo && sum(s)<=hi) acc <- rbind(acc, c(m, or_of(dd$cl,s)))
  if (!is.null(acc) && sum(acc[,1]==m) >= 1000) break }
for (m in mults) { a <- acc[acc[,1]==m,2,drop=TRUE]; if (length(a)) cat(sprintf("rate x%-3d accepted %4d | null OR median %.2f p95 %.2f | p = %.4f\n", m, length(a), median(a), quantile(a,.95), (sum(a>=obs)+1)/(length(a)+1))) }
a <- acc[,2]; cat(sprintf("Pooled conditioned null: %d replicates, median OR %.2f, p95 %.2f, p = %.4f\n", length(a), median(a), quantile(a,.95), (sum(a>=obs)+1)/(length(a)+1)))
## 3. Within-FSM opsin ~ genome size, phylogenetically corrected
f <- d[d$cl=="FSM" & !is.na(d$opsin) & d$tip %in% tr$tip.label,]
ft <- drop.tip(tr, setdiff(tr$tip.label, f$tip)); rownames(f) <- f$tip; f <- f[ft$tip.label,]
f$comp <- as.numeric(f$completeness_used_pct)
g1 <- glm(opsin ~ gs, data=f, family=binomial); cat("\nWithin FSM (n=",nrow(f),") plain logistic: beta ", round(coef(g1)[2],3), " p ", signif(summary(g1)$coef[2,4],3), "\n")
p1 <- phyloglm(opsin ~ gs, data=f, phy=ft, method="logistic_MPLE", btol=30)
print(summary(p1)$coefficients)
p2 <- phyloglm(opsin ~ gs + comp, data=f, phy=ft, method="logistic_MPLE", btol=30)
print(summary(p2)$coefficients); cat("alpha", p1$alpha, p2$alpha, "\n")
cat("Pagel lambda of opsin within FSM: "); print(phylosig(ft, setNames(f$opsin, f$tip), method="lambda"))
sink()
