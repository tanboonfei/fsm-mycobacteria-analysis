#!/bin/bash
# usage: relax_savenull.sh <alignment> <labelled tree> <out json> <threads>   (also writes <out>.fit.bf and <out>.fit.bf.null.bf)
A=$1; T=$2; O=$3; N=$4
OMP_NUM_THREADS=$N /usr/lib/hyphy/bin/hyphy-avx LIBPATH=/usr/share/hyphy CPU=$N /path/to/project/analyses/extended/D_relax/scripts/RELAX_savenull.bf --alignment $A --tree $T \
   --test Test --reference Reference --code Universal --models Minimal --srv No --output $O --save-fit ${O%.json}.fit.bf > ${O%.json}.log 2>&1
