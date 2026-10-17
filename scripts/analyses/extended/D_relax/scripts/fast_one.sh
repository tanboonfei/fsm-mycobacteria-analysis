#!/bin/bash
# usage: fast_one.sh <alignment|REAL> <tag>
D=/path/to/project/analyses/extended/D_relax
mkdir -p $D/boot/fast
a=$1; [ "$a" != "REAL" ] && [ ! -s "$a" ] && { echo "missing $a"; exit 1; }
printf "%s\n%s\n" $D/boot/real/FSM_tips_stem.fit.bf.null.bf $a | OMP_NUM_THREADS=1 /usr/lib/hyphy/bin/hyphy-avx LIBPATH=/usr/share/hyphy CPU=1 $D/scripts/fastrelax.bf > $D/boot/fast/$2.log 2>&1
