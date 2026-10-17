suppressPackageStartupMessages({library(ape); library(castor)})
B <- "/path/to/project/"; OUT <- paste0(B,"revision_checks/out/")
mat <- read.delim(paste0(OUT,"m1_presence_matrix.tsv"), row.names=1, check.names=FALSE)
meta<- read.delim(paste0(OUT,"m1_genome_meta.tsv"), stringsAsFactors=FALSE)
fates<-read.delim(paste0(OUT,"m1_family_fates.tsv"), stringsAsFactors=FALSE)
tr  <- read.tree(paste0(B,"iqtree/myco_bac120.rooted.treefile"))
run <- function(minc, tag){
  keep <- meta$genome[meta$tip %in% tr$tip.label & meta$completeness >= minc]
  md <- meta[match(keep,meta$genome),]; m <- mat[keep,,drop=FALSE]; rownames(m) <- md$tip
  tt <- drop.tip(tr, setdiff(tr$tip.label, md$tip)); m <- m[tt$tip.label,,drop=FALSE]
  md <- md[match(tt$tip.label, md$tip),]; ntip <- Ntip(tt)
  crown <- getMRCA(tt, md$tip[md$clade=="FSM"]); parent <- tt$edge[tt$edge[,2]==crown,1]
  P <- matrix(NA_real_, ncol(m), 2)
  for (j in seq_len(ncol(m))) { x <- m[,j]+1L
    if (length(unique(x))==1L) { P[j,] <- as.numeric(unique(x)==2L); next }
    f <- tryCatch(asr_mk_model(tt,x,Nstates=2,rate_model="ER",include_ancestral_likelihoods=TRUE,reroot=TRUE), error=function(e) NULL)
    if (!is.null(f) && f$success) P[j,] <- c(f$ancestral_likelihoods[crown-ntip,2], f$ancestral_likelihoods[parent-ntip,2]) }
  lost53 <- colnames(m) %in% fates$family[fates$fate=="lost_cladewide"]
  for (th in list(c(.95,.05), c(.90,.10))) {
    s <- P[,2]>=th[1] & P[,1]<=th[2]
    cat(sprintf("%s (%d tips: FSM %d, RGM %d, SGM %d) threshold %.2f/%.2f: stem losses %d (of the 53 clade-wide: %d); expected (sum p_prec*(1-p_crown)) %.1f\n",
      tag, ntip, sum(md$clade=="FSM"), sum(md$clade=="RGM"), sum(md$clade=="SGM"), th[1], th[2], sum(s,na.rm=TRUE), sum(s & lost53, na.rm=TRUE), sum(P[,2]*(1-P[,1]), na.rm=TRUE)))
  }
}
sink("asr178_out.txt", split=TRUE)
run(0, "243-genome panel (as published)")
run(95, "178-genome panel (>=95% complete)")
sink()
