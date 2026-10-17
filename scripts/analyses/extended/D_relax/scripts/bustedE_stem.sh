#!/bin/bash
# usage: bustedE_stem.sh <tree with {Test} on the stem> <out json> <threads> <alignment>
OMP_NUM_THREADS=$3 /usr/lib/hyphy/bin/hyphy-avx LIBPATH=/usr/share/hyphy CPU=$3 busted --alignment $4 --tree $1 --branches Test --code Universal --srv No --error-sink Yes --output $2 > ${2%.json}.log 2>&1
