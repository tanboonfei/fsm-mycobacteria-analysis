#!/bin/bash
while pgrep -f chain2.sh >/dev/null; do sleep 30; done
for o in $(cat og_keep.txt); do [ -s trim/$o.trim.faa ] || ./align.sh $o; done
python3 concat.py > concat.log 2>&1
echo ALN_DONE
