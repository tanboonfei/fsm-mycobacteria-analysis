#!/bin/bash
while pgrep -f chain1.sh >/dev/null; do sleep 30; done
# scOU with terminal edges floored at 0.002 tree heights (ultrametric re-imposed by extending terminals then rescaling)
Rscript run_pem.R 824 log myco_bac120.chronos_l1.nwk 1 scOU 0.002 _floor > pem_824_log_floor.log 2>&1
Rscript pem_post.R pem_824_log_floor.rds > pem_824_log_floor_post.log 2>&1
# BM shift model on the original phylogram
Rscript run_pem.R 824 log myco_bac120.rooted.treefile 1 BM 0 _BMphylogram > pem_824_log_BM.log 2>&1
Rscript pem_post.R pem_824_log_BMphylogram.rds > pem_824_log_BM_post.log 2>&1
