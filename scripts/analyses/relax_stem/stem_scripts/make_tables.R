## Replacement rows for Tables S6f / S6g from joint_stem_summary.tsv and per-family results.
DST <- "/path/to/project/analyses/relax_stem/stem/out/"
TB <- "/path/to/project/analyses/relax_stem/tables/"
s <- read.delim(paste0(DST, "joint_stem_summary.tsv"), check.names = FALSE)
d <- read.delim(gzfile(paste0(DST, "joint_stem_per_family.tsv.gz")))
sat_fun <- function(e) abs(e$castor_p_parent - 0.5) < 1e-3 & abs(e$castor_p_crown - 0.5) < 1e-3
## S6f: pooled-rate models, both panels
f <- s[s$model != "per-family ER (castor)", ]
S6f <- data.frame(Panel = ifelse(f$panel == "243", "243 genomes", "178 genomes (>=95% complete)"), Model = f$model,
  `Gain rate` = f$gain, `Loss rate` = f$loss,
  `Expected stem losses, joint P(present at parent, absent at crown), all 1,245` = f$joint_all,
  `Of which among the 53` = f$joint_53,
  `Families with joint P >= 0.95 (of which in 53)` = sprintf("%d (%d)", f$joint_ge95, f$joint_ge95_in53),
  `Families with joint P >= 0.90 (of which in 53)` = sprintf("%d (%d)", f$joint_ge90, f$joint_ge90_in53),
  `Expected loss events on stem (all / 53)` = sprintf("%.1f / %.1f", f$Eloss_stem_all, f$Eloss_stem_53),
  `Previously reported product sum (castor marginals)` = f$old_product_all,
  `Previously reported confident 0.95/0.05; 0.90/0.10` = sprintf("%d; %d", f$old_conf95, f$old_conf90), check.names = FALSE)
write.table(S6f, paste0(TB, "S6f_replacement_rows.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
## S6g: per-family ER
g <- s[s$model == "per-family ER (castor)", ]
rows <- list()
for (i in seq_len(nrow(g))) {
  e <- d[d$panel == g$panel[i] & d$model == "per-family ER (castor)", ]; sat <- sat_fun(e)
  rows[[i]] <- data.frame(Panel = ifelse(g$panel[i] == "243", "243-genome panel (published)", "178 genomes >=95% complete"),
    `Tips (FSM / RGM / SGM)` = ifelse(g$panel[i] == "243", "96 / 97 / 50", "59 / 83 / 36"),
    `Confident stem losses, marginal 0.95/0.05 (in 53) [unchanged]` = sprintf("%d (%d)", g$old_conf95[i], g$old_conf95_in53[i]),
    `Joint P >= 0.95 (in 53)` = sprintf("%d (%d)", g$joint_ge95[i], g$joint_ge95_in53[i]),
    `Joint P >= 0.90 (in 53)` = sprintf("%d (%d)", g$joint_ge90[i], g$joint_ge90_in53[i]),
    `Expected stem losses, sum of joint P, all 1,245` = g$joint_all[i],
    `Sum of joint P, excluding saturated fits` = sprintf("%.1f (n = %d excluded)", sum(e$joint_p10[!sat]), sum(sat)),
    `Sum of joint P, 53 clade-wide families` = g$joint_53[i],
    `Expected loss events on stem, 53 families (stochastic-mapping expectation)` = g$Eloss_stem_53[i] ,
    `Previously reported expected (product of marginals)` = g$old_product_all[i], check.names = FALSE)
}
S6g <- do.call(rbind, rows)
## 53-family E[loss events] excluding saturated (the 178-panel cwsA fit is saturated)
for (i in seq_len(nrow(g))) { e <- d[d$panel == g$panel[i] & d$model == "per-family ER (castor)" & d$fate == "lost_cladewide", ]; sat <- sat_fun(e)
  S6g[i, "Expected loss events on stem, 53 families (stochastic-mapping expectation)"] <- sprintf("%.1f%s", sum(e$E_loss_events_stem[!sat]), if (any(sat)) sprintf(" (excl. %d saturated)", sum(sat)) else "") }
write.table(S6g, paste0(TB, "S6g_replacement_rows.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
## per-family joint P for the 53 (243 panel, per-family ER and pooled ER), for an optional S6h column
e <- d[d$panel == "243" & d$fate == "lost_cladewide" & d$model %in% c("per-family ER (castor)", "pooled ER (gain = loss)", "gain = loss/10", "gain = loss/100 (Dollo-like)"), ]
w <- reshape(e[, c("family", "model", "joint_p10")], idvar = "family", timevar = "model", direction = "wide")
e178 <- d[d$panel == "178" & d$fate == "lost_cladewide" & d$model == "per-family ER (castor)", ]
w$`joint_p10.178 per-family ER` <- e178$joint_p10[match(w$family, e178$family)]
old <- d[d$panel == "243" & d$model == "per-family ER (castor)", ]
w$old_confident_33 <- old$old_conf95[match(w$family, old$family)]
w$prev_FSM_hc_pct <- old$prev_FSM_hc[match(w$family, old$family)]
write.table(w[order(-w$`joint_p10.per-family ER (castor)`), ], paste0(TB, "S6h_joint_probability_53.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(S6f); print(S6g)
