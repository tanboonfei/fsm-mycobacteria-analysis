for i in 0 1 2 3 4 5; do
 hmmsearch --cpu 1 -T 15 --noali -o /dev/null --tblout /path/to/scratch/kofam/part$i.tbl /path/to/scratch/kofam/part$i.hmm /path/to/scratch/kofam/sub.faa &
done
wait
echo done
