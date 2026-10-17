#!/bin/bash
# usage: relax_seed.sh <alignment> <labelled tree> <out json> <threads> <seed> <starting-points>
A=$1; T=$2; O=$3; N=$4; SEED=$5; SP=$6
S=/path/to/project/analyses/relax_stem/relax/scripts
/usr/bin/time -v -o ${O%.json}.time env OMP_NUM_THREADS=$N /usr/lib/hyphy/bin/hyphy-avx LIBPATH=/usr/share/hyphy CPU=$N ENV="multistart_seed=$SEED" $S/RELAX_seed.bf --alignment $A --tree $T \
   --test Test --reference Reference --code Universal --models Minimal --srv No --starting-points $SP --output $O --save-fit ${O%.json}.fit.bf > ${O%.json}.log 2>&1
