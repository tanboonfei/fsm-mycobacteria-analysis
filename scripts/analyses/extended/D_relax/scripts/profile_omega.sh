#!/bin/bash
# usage: profile_omega.sh <FSMstem|LEPstem> <tree> <threads>  -> full profile likelihood of the test-branch omega (grid in ci/<name>_prof.codon.fa.grid)
D=/path/to/project/analyses/extended/D_relax
printf "1\n$D/ci/$1_prof.codon.fa\n$2\n5\n" | OMP_NUM_THREADS=$3 /usr/lib/hyphy/bin/hyphy-avx LIBPATH=/usr/share/hyphy CPU=$3 $D/scripts/profile_omega.bf > $D/ci/profile_$1.log 2>&1
