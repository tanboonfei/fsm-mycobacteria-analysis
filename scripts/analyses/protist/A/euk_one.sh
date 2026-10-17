#!/bin/bash
# usage: euk_one.sh RUN URL NREADS OUTDIR  (streams reads; no temp files)
run=$1; url=$2; n=${3:-1000000}; OUT=${4:-k2out4}
D=/path/to/project/analyses/protist/A
mkdir -p $D/${OUT}; [ -s $D/${OUT}/$run.nreads ] && exit 0
source ~/miniconda3/etc/profile.d/conda.sh; conda activate kraken2_env
curl -s --retry 3 "$url" | zcat 2>/dev/null | head -n $((n*4)) | kraken2 --db $D/k2silva --threads 2 --confidence 0.1 --report $D/${OUT}/$run.report --output - /dev/stdin 2> $D/${OUT}/$run.log
nr=$(grep -oE '^[0-9]+ sequences' $D/${OUT}/$run.log | cut -d' ' -f1)
[ -n "$nr" ] && [ "$nr" -gt 1000 ] && echo -e "$run\t$nr" > $D/${OUT}/$run.nreads
