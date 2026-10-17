#!/usr/bin/env Rscript
# Draw the rooted opsin ML tree with mycobacterial tips coloured by clade and characterised references labelled
suppressPackageStartupMessages(library(ape))
O <- "/path/to/project/analyses/extended/B_opsin"
r <- ladderize(read.tree(file.path(O,"tree/opsin_ml.rooted_HeR.nwk")))
u <- read.tree(file.path(O,"tree/opsin_ml.treefile"))
A <- read.delim(file.path(O,"tree/tip_annotation.tsv"), stringsAsFactors=FALSE); rownames(A) <- A$tip
Lg <- read.delim(file.path(O,"screen/opsin_loci_lineages.tsv"), stringsAsFactors=FALSE); lin <- setNames(Lg$lineage, Lg$locus)
t <- r$tip.label; set <- A[t,"set"]; grp <- A[t,"group"]
col <- ifelse(set=="mycobacterial", ifelse(grp=="FSM", ifelse(lin[t]=="B","#b2182b","#2166ac"), ifelse(grp=="RGM","#e08214","#1b7837")),
        ifelse(set=="characterised","black", ifelse(A[t,"phylum"]=="Actinomycetota","#9970ab","grey60")))
lab <- ifelse(set=="characterised", A[t,"label"], ifelse(set=="mycobacterial", paste0(grp," ", t), paste0(A[t,"organism"]," (",A[t,"order"],")")))
pdf(file.path(O,"tree/opsin_ml_tree.pdf"), width=14, height=60)
plot(r, tip.color=col, cex=0.32, label.offset=0.02, show.tip.label=FALSE, no.margin=TRUE, edge.width=0.5)
tiplabels(lab, frame="none", adj=0, cex=0.32, col=col, offset=0.03)
legend("topleft", legend=c("FSM, lineage A","FSM, lineage B","RGM","SGM","characterised reference","other Actinomycetota (UniProt)","other UniProt"),
       text.col=c("#2166ac","#b2182b","#e08214","#1b7837","black","#9970ab","grey60"), bty="n", cex=1)
add.scale.bar(cex=0.8)
dev.off()
