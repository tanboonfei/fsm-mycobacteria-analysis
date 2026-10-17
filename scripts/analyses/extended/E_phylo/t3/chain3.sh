#!/bin/bash
for cfg in "BM rooted.treefile 1 1 _BM_derep_me" "scOU chronos_l1.nwk 1 1 _scOU_derep_me" "BM rooted.treefile 0 1 _BM_all_me" "scOU chronos_l1.nwk 0 1 _scOU_all_me"; do
 set -- $cfg
 Rscript run_pem.R 824 log myco_bac120.$2 1 $1 $3 $4 $5 > pem_824_log$5.log 2>&1
 Rscript pem_post.R pem_824_log$5.rds > pem_824_log${5}_post.log 2>&1
done
