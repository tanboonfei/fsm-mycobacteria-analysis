#!/bin/bash
# usage: warm.sh <saved fit> <start file|NONE> <null|alt> <out prefix> <condfirst 0|1> <threads>
S=/path/to/project/analyses/relax_stem/relax/scripts
printf '%s\n%s\n%s\n%s\n%s\n' "$1" "$2" "$3" "$4" "$5" | /usr/bin/time -v -o $4.time env OMP_NUM_THREADS=$6 /usr/lib/hyphy/bin/hyphy-avx LIBPATH=/usr/share/hyphy CPU=$6 $S/warm.bf > $4.log 2>&1
