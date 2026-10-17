suppressPackageStartupMessages({library(ape); library(nlme); library(phylolm)})
tr <- read.tree("inputs/myco_bac120.rooted.treefile")
d  <- read.delim("inputs/tsv/S1_genome_metadata.tsv",
                 stringsAsFactors=FALSE, na.strings=c("NA","","NaN","None"))
d$cl <- ifelse(grepl("^FSM",d$clade),"FSM", ifelse(grepl("^RGM",d$clade),"RGM",
        ifelse(d$clade=="SGM","SGM",NA)))
d <- d[!is.na(d$cl) & d$tip %in% tr$tip.label,]
tt <- drop.tip(tr, setdiff(tr$tip.label,d$tip)); rownames(d)<-d$tip; d<-d[tt$tip.label,]
d$FSM <- as.integer(d$cl=="FSM")
chk <- function(col,lab){
  dd <- d[!is.na(d[[col]]),]
  t2 <- drop.tip(tt, setdiff(tt$tip.label,dd$tip)); dd <- dd[t2$tip.label,]
  dd$y <- as.numeric(dd[[col]])
  g <- tryCatch(gls(y ~ FSM, data=dd, correlation=corPagel(0.9,phy=t2,form=~tip,fixed=FALSE),
                    method="ML"), error=function(e) NULL)
  if (is.null(g)) { cat(sprintf("%-22s gls FAILED\n",lab)); return(invisible()) }
  s <- summary(g)$tTable
  cat(sprintf("%-22s gls/corPagel(ML)  coef %+8.3f  lambda %.3f  p %.4g\n",
              lab, s["FSM","Value"], as.numeric(coef(g$modelStruct$corStruct, unconstrained=FALSE)), s["FSM","p-value"]))
}
chk("genome_size_Mb","genome size (Mb)")
chk("coding_density_pct_uniform_prodigal","coding density uniform")
chk("coding_density_pct","coding density deposited")
chk("acidic_residues_pct","acidic residues (%)")
chk("GC_pct","GC (%)")
