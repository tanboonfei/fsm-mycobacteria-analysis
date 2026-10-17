#!/bin/bash
# usage: profile.sh <saved fit> <K> <out prefix> <threads>
S=/path/to/project/analyses/relax_stem/relax/scripts
printf '%s\n%s\n%s\n' "$1" "$2" "$3" | /usr/bin/time -v -o $3.time env OMP_NUM_THREADS=$4 /usr/lib/hyphy/bin/hyphy-avx LIBPATH=/usr/share/hyphy CPU=$4 $S/profile.bf > $3.log 2>&1
