#!/bin/bash
# usage: relax_generic.sh <alignment> <labelled tree> <out json> <threads> [extra args...]
A=$1; T=$2; O=$3; N=$4; shift 4
OMP_NUM_THREADS=$N /usr/lib/hyphy/bin/hyphy-avx LIBPATH=/usr/share/hyphy CPU=$N relax --alignment $A --tree $T \
   --test Test --reference Reference --code Universal --models Minimal --srv No --output $O "$@" > ${O%.json}.log 2>&1
