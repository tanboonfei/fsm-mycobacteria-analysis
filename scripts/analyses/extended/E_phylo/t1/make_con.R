suppressMessages(library(ape))
s1 <- read.delim("../../../../Publication/scripts/inputs/tsv/S1_genome_metadata.tsv", check.names=FALSE)
cl <- setNames(as.character(s1$clade), s1$tip); sp <- setNames(as.character(s1$species), s1$tip)
tr <- read.tree("r1_lg.treefile"); og <- tr$tip.label[grepl("Outgroup", cl[tr$tip.label])]; tr <- root(tr, og, resolve.root=TRUE)
tr$node.label <- NULL
desc <- function(t,n) if (n<=Ntip(t)) t$tip.label[n] else extract.clade(t,n)$tip.label
fsm <- tr$tip.label[grepl("^FSM", cl[tr$tip.label])]
s42 <- readLines("sis_level1.txt"); pal <- readLines("sis_level2.txt"); l3 <- readLines("sis_level3.txt"); sgm <- readLines("sis_level4.txt")
m <- getMRCA(tr,s42); k <- tr$edge[tr$edge[,1]==m,2]; s42a <- desc(tr,k[1]); s42b <- desc(tr,k[2]); if (length(s42a)<length(s42b)) {t<-s42a;s42a<-s42b;s42b<-t}
m <- getMRCA(tr,l3); k <- tr$edge[tr$edge[,1]==m,2]; l3a <- desc(tr,k[1]); l3b <- desc(tr,k[2])
any4 <- tr$tip.label[sp[tr$tip.label] %in% c("Mycobacterium sp001494595","Mycobacterium anyangense","Mycobacterium sp011694515","Mycobacterium sp036023875")]
stopifnot(length(any4)==4, is.monophyletic(tr, any4))
rgm <- tr$tip.label[grepl("^RGM", cl[tr$tip.label])]
H <- list(H1_anyangense4=any4, H2_S42_rhodesiae30=s42a, H3_S42_sp020438785_12=s42b, H4_palauense4=pal, H5_S42plusPalauense46=c(s42,pal),
          H6_RGM_chubuense_phlei19=if ("Mycobacterium chubuense" %in% sp[l3a]) l3a else l3b, H7_RGM_fortuitum_neoaurum18=if ("Mycobacterium chubuense" %in% sp[l3a]) l3b else l3a,
          H8_SGM30=sgm, H9_allRGM83=rgm)
fsub <- extract.clade(tr, getMRCA(tr,fsm)); rest <- drop.tip(tr, fsm)
info <- data.frame()
for (h in names(H)) {
  X <- H[[h]]; stopifnot(is.monophyletic(tr, X) || h %in% c("H5_S42plusPalauense46","H9_allRGM83"))
  others <- setdiff(tr$tip.label, c(fsm,X))
  con <- sprintf("(((%s),(%s)),%s);", paste(fsm,collapse=","), paste(X,collapse=","), paste(others,collapse=","))
  writeLines(con, paste0("con_",h,".nwk"))
  n <- if (length(X)==1) which(rest$tip.label==X) else getMRCA(rest, X)
  rs <- write.tree(rest); xs <- sub(";$","",write.tree(extract.clade(rest,n))); fs <- sub(";$","",write.tree(fsub))
  stopifnot(grepl(xs, rs, fixed=TRUE)); st <- read.tree(text=sub(xs, paste0("(",xs,":0.01,",fs,":0.01)"), rs, fixed=TRUE))
  stopifnot(is.monophyletic(st, c(fsm,X)), is.monophyletic(st, fsm))
  write.tree(unroot(st), paste0("start_",h,".nwk"))
  info <- rbind(info, data.frame(hypothesis=h, n_X=length(X), X_species_examples=paste(head(unique(sp[X]),4),collapse="; "), X_monophyletic_in_ML=is.monophyletic(tr,X)))
}
write.table(info, "hypotheses.tsv", sep="\t", quote=FALSE, row.names=FALSE); print(info)
