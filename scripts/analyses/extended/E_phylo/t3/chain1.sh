#!/bin/bash
Rscript run_pem.R 824 log myco_bac120.chronos_l1.nwk 1 > pem_824_log.log 2>&1
Rscript run_l1ou.R 824 raw pBIC myco_bac120.chronos_l1.nwk > l1ou_824_raw_pBIC_chronos.log 2>&1
Rscript run_l1ou.R 525 log pBIC myco_bac120.chronos_l1.nwk > l1ou_525_log_pBIC_chronos.log 2>&1
Rscript run_l1ou.R 824 log BIC myco_bac120.chronos_l1.nwk > l1ou_824_log_BIC_chronos.log 2>&1
Rscript run_l1ou.R 824 log pBIC myco_bac120.extend.nwk > l1ou_824_log_pBIC_extend.log 2>&1
