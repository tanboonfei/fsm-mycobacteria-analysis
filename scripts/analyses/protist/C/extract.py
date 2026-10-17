import sys
seen=set(); out=open('all862.faa','w'); hdr=None; seq=[]
def flush():
    if hdr and seq:
        k=hdr.split()[0]
        if k not in seen:
            seen.add(k); out.write('>'+k+'\n'+''.join(seq)+'\n')
for line in open('../../../mmclust862/c862_all_seqs.fasta'):
    if line.startswith('>'):
        flush(); hdr=line[1:].strip(); seq=[]
    else: seq.append(line.strip())
flush(); print(len(seen))
