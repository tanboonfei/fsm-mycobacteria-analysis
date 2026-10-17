source("common.R"); suppressPackageStartupMessages(library(l1ou))
fs <- c("l1ou_824_log_pBIC_chronos_l1.rds","l1ou_824_log_BIC_chronos_l1.rds","l1ou_824_raw_pBIC_chronos_l1.rds","l1ou_525_log_pBIC_chronos_l1.rds","l1ou_824_log_pBIC_extend.rds","l1ou_824_log_pBIC_chronos_l1_derep.rds")
all <- data.frame(); reg <- data.frame()
for (f in fs) {
  x <- readRDS(f); fit <- x$fit; tt <- x$tree; lg <- grepl("_log_", f); tr <- function(v) if (lg) exp(v) else v
  sh <- fit$shift.configuration
  for (i in seq_along(sh)) all <- rbind(all, cbind(run=sub("\\.rds$","",f), nShifts=length(sh), alpha_unit_height=signif(fit$alpha,4), describe_edge(tt, sh[i], fit$Y[,1]),
        shift_in_optimum=signif(fit$shift.values[i],4), shift_in_tip_mean=signif(fit$shift.means[i],4), tip_mean_ratio=if (lg) signif(exp(fit$shift.means[i]),3) else NA))
  op <- fit$optima; if (is.matrix(op)) op <- op[,1]; names(op) <- tt$tip.label
  mu <- fit$mu; if (is.matrix(mu)) mu <- mu[,1]; names(mu) <- tt$tip.label
  for (g in c("FSM","sister42","leprosy")) { tp <- intersect(switch(g, FSM=fsm_tips, sister42=s42, leprosy=lep_tips), tt$tip.label); if (!length(tp)) next
    reg <- rbind(reg, data.frame(run=sub("\\.rds$","",f), group=g, n=length(tp), expected_mean_Mb_median=round(median(tr(mu[tp])),3), optimum_Mb_modal=round(tr(as.numeric(names(which.max(table(round(op[tp],6)))))),3), observed_median_Mb=round(median(d[tp,"genome_size_Mb"]),3))) }
}
write.table(all, "l1ou_shifts_all.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(reg, "l1ou_regimes.tsv", sep="\t", quote=FALSE, row.names=FALSE)
print(all[all$label!="",c("run","nShifts","n_tips","label","shift_in_optimum","shift_in_tip_mean","tip_mean_ratio")]); print(reg)
