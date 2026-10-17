tblastn -query /path/to/scratch/tbn/q.faa -db /path/to/project/analyses/protist/C/tbn/fsm -evalue 1e-5 -max_target_seqs 20000 -num_threads 6 -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" > /path/to/scratch/tbn/hits_fsm.tsv
tblastn -query /path/to/scratch/tbn/q.faa -db /path/to/project/analyses/protist/C/tbn2/sis -evalue 1e-5 -max_target_seqs 20000 -num_threads 6 -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" > /path/to/scratch/tbn/hits_sis.tsv
cmsearch --cut_ga --cpu 6 --noali --tblout /path/to/scratch/tbn/tpp_fsm.tbl -o /dev/null /path/to/scratch/rfam_dl/RF00059.cm /path/to/project/analyses/protist/C/tbn/fsm.fna
cmsearch --cut_ga --cpu 6 --noali --tblout /path/to/scratch/tbn/tpp_sis.tbl -o /dev/null /path/to/scratch/rfam_dl/RF00059.cm /path/to/project/analyses/protist/C/tbn2/sis.fna
echo ALLDONE
