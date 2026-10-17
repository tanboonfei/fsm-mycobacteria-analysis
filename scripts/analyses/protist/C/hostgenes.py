import collections,re
from groups import load
G=load()
SETS={
'Cholesterol catabolism':['kstD','kshA','kshB','hsaA','hsaB','hsaC','hsaD','cyp125','cyp142','fadA5','ipdA','ipdB'],
'Mce1':['mce1A','mce1B','mce1C','mce1D','mce1E','mce1F'],'Mce2':['mce2A','mce2B','mce2C','mce2D'],'Mce3':['mce3A','mce3B','mce3C','mce3D','mce3E','mce3F'],'Mce4':['mce4A','mce4B','mce4C','mce4D','mce4E','mce4F'],
'ESX-1':['eccCa1','eccCb1','eccD1','eccB1','eccE1','mycP1','espA','espB','espC'],'ESX-3':['eccC3','eccD3','eccB3','eccE3','mycP3'],'ESX-4':['eccC4','eccD4','eccB4','mycP4','esxT','esxU'],'ESX-5':['eccC5','eccD5','eccB5','mycP5'],
'Mycobactin siderophore':['mbtA','mbtB','mbtC','mbtD','mbtE','mbtF','mbtG','mbtH','mbtI','mbtJ'],'IrtAB siderophore import':['irtA','irtB'],
'Oxidative-stress defence':['katG','katE','sodA','sodC','ahpC','ahpD','trxB'],'Phagosome-arrest phosphatases':['ptpA','ptpB','sapM'],
'PDIM/PGL lipids':['ppsA','ppsB','ppsC','ppsD','ppsE','mas','drrA','drrB','fadD26','fadD28','mmpL7','papA5'],
'Ergothioneine':['egtA','egtB','egtC','egtD','egtE'],'Cytochrome bd':['cydA','cydB','cydC','cydD'],
}
names=collections.defaultdict(set)
for l in open('../../../eggnog_ref_out/allfam.emapper.annotations'):
    if l[0]=='#':continue
    f=l.split('\t'); 
    if f[8] not in ('-',''): names[f[0]]=set(x.strip() for x in f[8].split(','))
gn=collections.defaultdict(set)
for l in open('../../../mmclust862/c862_cluster.tsv'):
    r,m=l.split()
    if r in names: gn[m.split('|')[0]].update(names[r])
allg=set(x for v in SETS.values() for x in v)
found=collections.Counter(x for s in gn.values() for x in s if x in allg)
print('names never seen:',sorted(allg-set(found)))
order=['FSM','RGM_sister','RGM_other','SGM_env','NTM_opportunistic','Leprosy','SGM_obligate']
o=open('hostgene_modules.tsv','w'); o.write('module\tgenes_found\t'+'\t'.join(order)+'\n')
for mod,genes in SETS.items():
    gs=[x for x in genes if x in found]
    row=[mod,','.join(gs)]
    for grp in order:
        sub=[g for g,v in G.items() if v['group']==grp and v['comp']>=90]
        frac=[sum(x in gn[g] for x in gs)/len(gs) for g in sub] if gs else []
        row.append(f"{100*sum(frac)/len(frac):.0f}" if frac else 'NA')
    o.write('\t'.join(row)+'\n')
o.close()
# per gene too
o=open('hostgene_pergene.tsv','w'); o.write('gene\t'+'\t'.join(order)+'\n')
for mod,genes in SETS.items():
    for x in genes:
        if x not in found: continue
        o.write(x+'\t'+'\t'.join(f"{100*sum(x in gn[g] for g in [g for g,v in G.items() if v['group']==grp and v['comp']>=90])/max(1,sum(1 for v in G.values() if v['group']==grp and v['comp']>=90)):.0f}" for grp in order)+'\n')
