import os, pandas as pd
GROUPS={'Eukaryota':'D','Bacteria':'D','Archaea':'D','Animalia':'K','Chloroplastida':'K','Fungi':'K','Amoebozoa':'K','Ciliophora':'P','Cercozoa':'P','Discoba':None,'Kinetoplastea':None,'Heterolobosea':None,'Chrysophyceae':None,'Cryptophyceae':None,'Choanoflagellida':None,'Dinoflagellata':None,'Diatomea':'P','Rhizaria':'K','Alveolata':'K','Stramenopiles':'K','Apicomplexa':None,'Centrohelida':None,'Bicosoecida':None,'Telonemia':None,'Mycobacterium':'G','Tubulinea':None,'Discosea':None,'Evosea':None,'Rickettsiales':'O','Chlamydiales':'O','Legionellales':'O','Nucletmycea':None,'Haptophyta':None,'Picozoa':None}
def parse(fn):
    d={k:0 for k in GROUPS}; seen=set()
    for line in open(fn):
        f=line.rstrip('\n').split('\t')
        if len(f)<6: continue
        name=f[5].strip()
        if name in GROUPS and name not in seen:
            d[name]=int(f[1]); seen.add(name)
    return d
def table(runs,outdir='k2out4'):
    rows=[]
    for r in runs:
        rep=f'{outdir}/{r}.report'; nr=f'{outdir}/{r}.nreads'
        if os.path.exists(rep) and os.path.exists(nr):
            d=parse(rep); d['run_accession']=r; d['nreads']=int(open(nr).read().split()[1]); rows.append(d)
    t=pd.DataFrame(rows)
    if len(t):
        t=t.rename(columns={c:'k2_'+c for c in ['Rickettsiales','Chlamydiales','Legionellales','Mycobacterium']})
        t['SSU_total']=t.Bacteria+t.Archaea+t.Eukaryota
        # protists = eukaryotes minus animals, plants (Chloroplastida), fungi
        t['protists']=t.Eukaryota-t.Animalia-t.Chloroplastida-t.Fungi
    return t
