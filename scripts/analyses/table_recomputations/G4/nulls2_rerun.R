suppressPackageStartupMessages({library(ape);library(phytools);library(phylolm)})
tr <- read.tree("/path/to/project/scripts/inputs/myco_bac120.rooted.treefile")
d  <- read.delim("/path/to/project/scripts/inputs/tsv/S1_genome_metadata.tsv", stringsAsFactors=FALSE, na.strings=c("NA","","NaN","None"))
d$cl <- ifelse(grepl("^FSM",d$clade),"FSM", ifelse(grepl("^RGM",d$clade),"RGM", ifelse(d$clade=="SGM","SGM",NA)))
d <- d[!is.na(d$cl) & !duplicated(d$tip),]; d$gs <- as.numeric(d$genome_size_Mb); d$opsin <- ifelse(is.na(d$rhodopsin), NA, as.integer(d$rhodopsin=="present"))
sink("rerun/nulls2_out_new.txt", split=TRUE)
## within-FSM opsin ~ genome size, phylogenetically corrected
f <- d[d$cl=="FSM" & !is.na(d$opsin) & d$tip %in% tr$tip.label,]
ft <- drop.tip(tr, setdiff(tr$tip.label, f$tip)); rownames(f) <- f$tip; f <- f[ft$tip.label,]; f$comp <- as.numeric(f$completeness_used_pct)
g1 <- glm(opsin ~ gs, data=f, family=binomial); cat("Within FSM n=",nrow(f)," plain logistic beta", round(coef(g1)[2],3), "p", signif(summary(g1)$coef[2,4],3), "\n")
p1 <- phyloglm(opsin ~ gs, data=f, phy=ft, method="logistic_MPLE", btol=30); print(summary(p1)$coefficients); cat("alpha", p1$alpha, "\n")
p2 <- phyloglm(opsin ~ gs + comp, data=f, phy=ft, method="logistic_MPLE", btol=30); print(summary(p2)$coefficients)
s <- f[f$origin=="SPIRE MAG",]; st <- drop.tip(ft, setdiff(ft$tip.label, s$tip)); s <- s[st$tip.label,]
p3 <- phyloglm(opsin ~ gs, data=s, phy=st, method="logistic_MPLE", btol=30); cat("SPIRE MAGs only n=", nrow(s), "\n"); print(summary(p3)$coefficients)
## prevalence-conditioned opsin null under fitted ARD model
set.seed(20260901)
or_of <- function(cl, y){ tab <- table(factor(cl,levels=c("RGM","FSM")), factor(y,levels=c(0,1)))+0.5; (tab[2,2]*tab[1,1])/(tab[2,1]*tab[1,2]) }
dd <- d[d$cl %in% c("FSM","RGM") & !is.na(d$opsin) & d$tip %in% tr$tip.label,]
tt <- drop.tip(tr, setdiff(tr$tip.label, dd$tip)); rownames(dd) <- dd$tip; dd <- dd[tt$tip.label,]
y <- setNames(as.character(dd$opsin), dd$tip); npos <- sum(dd$opsin); obs <- or_of(dd$cl, dd$opsin)
fit <- fitMk(tt, y, model="ARD"); Q <- as.Qmatrix(fit); print(Q)
pi1 <- Q["0","1"]/(Q["0","1"]+Q["1","0"]); cat("stationary P(opsin) =", round(pi1,3), " observed prevalence =", round(npos/length(y),3), "\n")
lo <- floor(npos*0.9); hi <- ceiling(npos*1.1); acc <- c(); tried <- 0
while (length(acc) < 2000 && tried < 60000) { tried <- tried+1
  r <- sim.Mk(tt, Q, anc=sample(c("0","1"),1,prob=c(1-pi1,pi1)))
  s <- as.integer(as.character(r[tt$tip.label])); if (sum(s)>=lo && sum(s)<=hi) acc <- c(acc, or_of(dd$cl, s)) }
cat(sprintf("Conditioned ARD null: tried %d, accepted %d (positives %d-%d); null OR median %.2f, p95 %.2f, max %.2f; observed %.2f; p = %.4f\n", tried, length(acc), lo, hi, median(acc), quantile(acc,.95), max(acc), obs, (sum(acc>=obs)+1)/(length(acc)+1)))
sink()
## opsin phyloglm ~ clade (FSM vs RGM) with current S1
sink("rerun/opsin_phyloglm.txt")
dd$clade <- factor(dd$cl, levels=c("RGM","FSM")); dd$y <- dd$opsin
print(table(dd$clade, dd$y)); ft <- fisher.test(table(dd$clade, dd$y)); print(ft)
f2 <- phyloglm(y ~ clade, data=dd, phy=tt, method="logistic_MPLE", btol=50); print(summary(f2)$coefficients); cat("alpha", f2$alpha, "\n")
sink()
