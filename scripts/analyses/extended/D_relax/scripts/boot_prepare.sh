#!/bin/bash
# simulate N replicates under the fitted RELAX null (K:=1) for the FSM terminals+stem labelling; impose real missing-data pattern
D=/path/to/project/analyses/extended/D_relax; N=$1
mkdir -p $D/boot/sim; cd $D/boot/sim
printf "%s\n%s\n%s\n" $D/boot/real/FSM_tips_stem.fit.bf.null.bf $D/boot/sim/raw $N | /usr/lib/hyphy/bin/hyphy-avx LIBPATH=/usr/share/hyphy CPU=1 $D/scripts/simulate.bf > $D/boot/sim/simulate.log 2>&1
for i in $(seq 1 $N); do python3 $D/scripts/postsim.py $D/boot/sim/raw.$i.fa $D/inputs/three_way.codon.fa $D/boot/sim/rep$i.fa && rm -f $D/boot/sim/raw.$i.fa; done
