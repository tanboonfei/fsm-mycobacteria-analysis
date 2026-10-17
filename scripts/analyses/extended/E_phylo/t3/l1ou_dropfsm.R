# information-criterion cost of removing the FSM-stem shift from each l1ou configuration (other shifts fixed, alpha re-estimated)
source("common.R"); suppressPackageStartupMessages(library(l1ou))
fs <- c("l1ou_824_log_pBIC_chronos_l1.rds","l1ou_824_raw_pBIC_chronos_l1.rds","l1ou_525_log_pBIC_chronos_l1.rds","l1ou_824_log_pBIC_chronos_l1_derep.rds")
out <- data.frame()
for (f in fs) { if (!file.exists(f)) next
  x <- readRDS(f); tt <- x$tree; Y <- x$Y; sh <- x$fit$shift.configuration
  fe <- which(tt$edge[,2]==getMRCA(tt, intersect(fsm_tips, tt$tip.label)))
  lp <- intersect(lep_tips, tt$tip.label)
  for (crit in c("BIC","AICc","pBIC")) {
    full <- fit_OU(tt, Y, sh, criterion=crit); red <- fit_OU(tt, Y, setdiff(sh, fe), criterion=crit)
    out <- rbind(out, data.frame(run=sub("\\.rds$","",f), criterion=crit, FSM_stem_in_config=fe %in% sh, score_full=round(full$score,2), score_without_FSM_shift=round(red$score,2), delta=round(red$score-full$score,2),
                 logLik_full=round(full$logLik,2), logLik_without=round(red$logLik,2), alpha_full=signif(full$alpha,3), alpha_without=signif(red$alpha,3)))
  }
}
print(out); write.table(out, "l1ou_drop_fsm_shift.tsv", sep="\t", quote=FALSE, row.names=FALSE)
