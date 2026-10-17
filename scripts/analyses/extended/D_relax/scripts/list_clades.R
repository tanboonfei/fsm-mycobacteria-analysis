suppressMessages(library(ape)); library(jsonlite)
B <- "/path/to/project/"
x <- readRDS("clades.rds"); pt <- x$pt
s1 <- read.delim(paste0(B,"Publication/scripts/inputs/tsv/S1_genome_metadata.tsv"), check.names=FALSE, na.strings=c("NA","","NaN","None"))
s1$id <- sub("\\.fa$","",s1$tip)
out <- list(); i <- 0
for (c in x$maxi) { ss <- s1[match(c$tips, s1$id),]
  ss <- ss[ss$genome_size_Mb >= 4.0,]
  if (nrow(ss) < 3) next
  i <- i+1; nm <- sprintf("CC%02d_%s", i, ifelse(grepl("^RGM",ss$clade[1]),"RGM","SGM"))
  out[[nm]] <- list(tips=ss$id, species=unique(ss$species), size_Mb=ss$genome_size_Mb, compl=ss$completeness_used_pct, dropped=setdiff(c$tips, ss$id)) }
for (n in names(out)) cat(n, length(out[[n]]$tips), "tips;", paste(out[[n]]$species, collapse="; "), "; dropped:", paste(out[[n]]$dropped, collapse=","), "\n")
write_json(out, "control_clades.json", pretty=TRUE, auto_unbox=FALSE)
