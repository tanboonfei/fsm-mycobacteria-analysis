import pandas as pd,re,pickle
d=pd.read_pickle('s1.pkl'); anno=pickle.load(open('anno.pkl','rb'))
sp=d[d.origin=='SPIRE MAG'].copy(); sp['id']=sp.tip.str.replace('.fa','',regex=False)
grp={'RGM':sp[sp.clade.astype(str).str.startswith('RGM')].id.tolist(),
     'FSM':sp[sp.clade.astype(str).str.startswith('FSM')].id.tolist(),
     'Leprosy clade':['GCF_000195855.1','GCF_000975265.2','GCF_003408705.1']}
R={ # gene symbol regex, product regex
 'ureC':(r'^ureC$', r'^urease( subunit alpha| alpha subunit)?$'),
 'narG':(r'^narG$', r'^(respiratory )?nitrate reductase( [0-9])?( subunit)? alpha( subunit| chain| apoprotein)?( apoprotein)?$'),
 'catA':(r'^catA$', r'^catechol 1,2-dioxygenase$'),
 'ectC':(r'^ectC$', r'^L-ectoine synthase$'),
 'kdpB':(r'^kdpB$', r'^(potassium|K\(\+\))-transporting ATPase (ATP-binding subunit|subunit B|subunit KdpB|B chain|KdpB subunit)$'),
}
def pos(t,gs,ps):
    m=(t.gene.str.fullmatch(gs.strip('^$')))|(t['product'].str.contains(ps,case=False,regex=True))
    return bool(m.any())
if __name__=='__main__':
    res={}
    for gname,ids in grp.items():
        ids2=[i for i in ids if i in anno]
        print(gname,len(ids),len(ids2),[i for i in ids if i not in anno])
        for k,(gs,ps) in R.items():
            p=[i for i in ids2 if pos(anno[i],gs,ps)]
            sym=[i for i in ids2 if anno[i].gene.str.match(gs).any()]
            print('  ',k,len(p),round(100*len(p)/len(ids2),1),'symbol-only',len(sym))
            res[(gname,k)]=p
    pickle.dump(res,open('res.pkl','wb'))
