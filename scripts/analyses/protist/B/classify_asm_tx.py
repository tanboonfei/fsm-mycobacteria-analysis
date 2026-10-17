"""Classify EukProt assembled-transcriptome contigs that hit FSM/control markers or 16S at >=90%.
16S: best FSM vs best non-FSM panel 16S (>=200 bp aligned); FSM if FSM>=98.7 and >= non-FSM+0.5.
Protein-coding: DIAMOND blastx vs 862-genome proteins; FSM if best FSM bitscore > best non-FSM and identity>=95
over >=50 aa."""
import json,collections
grp=json.load(open('panel_groups.json'))
fsm=set(x.replace('GB_','').replace('RS_','') for x in open('fsm_ids.txt').read().split())
fsm16=set(open('fsm_ids.txt').read().split())
def isF16(s): return grp.get(s,'').startswith('FSM') or s.split('_16S')[0] in fsm16
bad={'spire_mag_01108461','spire_mag_01753944','spire_mag_01960731','spire_mag_02102612'}
r16=collections.defaultdict(lambda:[0.0,0.0,'',0])
for l in open('asm_tx_vs_panel16S.tsv'):
    q,s,p,L,ql=l.split('\t'); p=float(p); L=int(L)
    if L<200 or s.split('_16S')[0] in bad: continue
    k=0 if isF16(s) else 1
    if p>r16[q][k]:
        r16[q][k]=p
        if k==1: r16[q][2]=s+':'+grp.get(s,'')
    r16[q][3]=max(r16[q][3],L)
rp=collections.defaultdict(lambda:[0.0,0.0,0.0,0.0,''])
for l in open('asm_tx_vs_myco.tsv'):
    q,s,p,L,ql,sl,e,b=l.split('\t'); p=float(p); b=float(b); g=s.split('|')[0]
    if int(L)<50: continue
    if g in fsm:
        if b>rp[q][1]: rp[q][0],rp[q][1]=p,b
    elif b>rp[q][3]: rp[q][2],rp[q][3],rp[q][4]=p,b,s
out=open('eukprot_asm_tx_classified.tsv','w')
out.write('transcript\tEukProt_ID\ttype\tbest_FSM_pid\tbest_nonFSM_pid\tbest_nonFSM\taln_len_or_bits\tcall\n')
cnt=collections.Counter()
for q in open('asm_hit_tx.txt').read().split():
    ep=q.split('|')[0]
    if q in r16 and max(r16[q][:2])>=95:
        F,N,h,L=r16[q]; t='16S'
        call='FSM' if F>=N+0.5 and F>=98.7 else ('nonFSM_myco' if N>=F else 'ambiguous')
        out.write(f'{q}\t{ep}\t16S\t{F}\t{N}\t{h}\t{L}\t{call}\n')
    elif q in rp and max(rp[q][0],rp[q][2])>=90:
        Fp,Fb,Np,Nb,h=rp[q]; t='protein'
        call='FSM' if (Fb>Nb and Fp>=95) else ('nonFSM_myco' if Nb>=Fb else 'ambiguous')
        out.write(f'{q}\t{ep}\tprotein\t{Fp}\t{Np}\t{h}\t{max(Fb,Nb)}\t{call}\n')
    else: continue
    cnt[(t,call)]+=1
print(cnt)
