#!/bin/bash
run=$1; url=$2; n=${3:-4000000}; OUT=${4:-k2out4}
D=/path/to/project/analyses/protist/A
for a in 1 2 3; do
  [ -s $D/$OUT/$run.nreads ] && exit 0
  rm -f $D/$OUT/$run.report
  $D/euk_one.sh $run $url $n $OUT
done
