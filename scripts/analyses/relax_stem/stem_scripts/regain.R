## 'None is regained' check for the 33 confident stem losses (243-genome panel, per-family ER, castor 0.95/0.05).
suppressPackageStartupMessages(library(ape))
DST <- "/path/to/project/analyses/relax_stem/stem/out/"
d <- read.delim(gzfile(paste0(DST, "joint_stem_per_family.tsv.gz")))
P <- readRDS(paste0(DST, "panel_243.rds")); P178 <- readRDS(paste0(DST, "panel_178.rds"))
fsm <- P$md$tip[P$md$clade == "FSM"]; fsm_hc <- P178$md$tip[P178$md$clade == "FSM"]
e <- d[d$panel == "243" & d$model == "per-family ER (castor)", ]
c33 <- e$family[e$old_conf95]
out <- data.frame(family = c33,
  n_FSM_present_243 = colSums(P$X[fsm, c33]), n_FSM_hc_present = colSums(P$X[fsm_hc, c33]), n_FSM_hc = length(fsm_hc),
  prev_FSM_hc_pct = e$prev_FSM_hc[match(c33, e$family)])
for (m in c("per-family ER (castor)", "pooled ER (gain = loss)", "gain = loss/10", "gain = loss/100 (Dollo-like)")) {
  z <- d[d$panel == "243" & d$model == m, ]; z <- z[match(c33, z$family), ]
  out[[paste0("joint_p10 | ", m)]] <- z$joint_p10
  out[[paste0("E_gains_in_FSM | ", m)]] <- z$E_gain_events_FSMcrown
  out[[paste0("E_losses_in_FSM | ", m)]] <- z$E_loss_events_FSMcrown
}
out <- out[order(-out$n_FSM_hc_present, out$family), ]
write.table(out, paste0(DST, "regain_33.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
det <- out$prev_FSM_hc_pct >= 1.5
cat(sprintf("33 confident families: detected in >=1 of %d high-completeness FSM MAGs: %d (prev 2-8%%: %d); detected in >=1 of %d FSM MAGs (243 panel): %d\n",
  length(fsm_hc), sum(out$n_FSM_hc_present > 0), sum(det), length(fsm), sum(out$n_FSM_present_243 > 0)))
print(table(out$n_FSM_hc_present))
for (m in c("per-family ER (castor)", "pooled ER (gain = loss)", "gain = loss/10", "gain = loss/100 (Dollo-like)")) {
  g <- out[[paste0("E_gains_in_FSM | ", m)]]
  cat(sprintf("%-30s expected gain events inside the FSM crown, summed over the 33: %.1f (over the %d detected families: %.1f); families with E[gains] >= 1: %d\n",
    m, sum(g), sum(out$n_FSM_hc_present > 0 | out$n_FSM_present_243 > 0), sum(g[out$n_FSM_hc_present > 0 | out$n_FSM_present_243 > 0]), sum(g >= 1)))
}
print(out[, 1:8], row.names = FALSE)
