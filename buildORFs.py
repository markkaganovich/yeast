import os
dir = './data/Alignments/ScerKwal/'
kwalalign = {}
for f in os.listdir(dir):
    file = open(dir+f)
    lines = file.readlines()
    file.close()
    f = f.split('.')[0]
    kwalalign[f] = {}
    for line in lines:    
        if line.startswith('>'):
            species = line.strip('\n').split('>')[1]
            seq = ''
        else:
            seq = seq + line.strip('\n')
        kwalalign[f][species] = seq




# build kwal alignment index
for o in orfs.keys():
    orf = orfs[o]
    kwalalignindex = {}
    if hasattr(orf, 'kwalalign') and 'Scer' in orf.kwalalign.keys() and 'Kwal' in orf.kwalalign.keys():
        dashes = 0
        i = 0
        while i < len(orf.seq):
            if orf.kwalalign['Scer'][i+dashes] == '-':
                dashes = dashes + 1
            else:
                kwalalignindex[i] = dashes + i
                i = i+1
        setattr(orf, 'kwalalignindex', kwalalignindex)

# build other alignments
otheralign = 'scasalign'
for o in orfs.keys():
    orf = orfs[o]
    otheralignindex = {}
    if hasattr(orf, otheralign) and 'Scer' in getattr(orf, otheralign).keys() and 'Scas' in getattr(orf, otheralign).keys():
        dashes = 0
        i = 0
        while i < len(orf.seq):
            if getattr(orf, otheralign)['Scer'][i+dashes] == '-':
                dashes = dashes + 1
            else:
                otheralignindex[i] = dashes + i
                i = i+1
        setattr(orf, otheralign+'index', otheralignindex)

# build other alignments




# phosphosite data sanity check
toolongfororf = []
notphospho = []
for o in orfs.keys():
    orf = orfs[o]
    if hasattr(orf, 'phosphosite'):
        l = map(lambda x: len(orf.seq) < x, orf.phosphosite)
        if sum(l) > 0:
            toolongfororf.append(orf.name)
            setattr(orf, 'phosphosite',[])
            continue
        t = filter(lambda x: orf.seq[x] in phosphorylatedaa, orf.phosphosite)
        if len(t)!=len(orf.phosphosite):
            notphospho.append(orf.name)
            setattr(orf, 'phosphosite',[])


# ortho and phosphosite conservation 
for o in orfs.keys():
    orf = orfs[o]
    if hasattr(orf, 'kwalalign') and 'Scer' in orf.kwalalign.keys() and 'Kwal' in orf.kwalalign.keys():
        (setattr(orf,'seqcons', 
                 sum(map(lambda x,y: x == y, orf.kwalalign['Scer'], orf.kwalalign['Kwal']))/
                 float(len(orf.seq))))
        if hasattr(orf, 'phosphosite') and len(orf.phosphosite) > 0:
            phosphoindex = map(lambda x: orf.kwalalignindex[x], orf.phosphosite)
            pcons = map(lambda x: orf.kwalalign['Scer'][x] == orf.kwalalign['Kwal'][x], phosphoindex)
            setattr(orf, 'phosphocons', sum(pcons)/float(len(orf.phosphosite)))

# syt conservation
for o in orfs.keys():
    orf = orfs[o]
    if hasattr(orf, 'kwalalign') and 'Scer' in orf.kwalalign.keys() and 'Kwal' in orf.kwalalign.keys():
        sytpos = ([i for i,x in enumerate(orf.seq) 
                if x == 'S' or x == 's' or x == 'T' or x == 't' or x == 'Y' or x == 'y'])
        sytposalign = map(lambda x: orf.kwalalignindex[x], sytpos)
        sytcons = map(lambda x: orf.kwalalign['Scer'][x] == orf.kwalalign['Kwal'][x], sytposalign)
        setattr(orf, 'sytcons', sum(sytcons)/float(len(sytpos)))



# add gene names to ORFs names
file = open('./data/go_slim_mapping.txt')
lines = file.readlines()
file.close()
for line in lines:
    if line.split('\t')[0] in orfs.keys():
        setattr(orfs[line.split('\t')[0]], 'gene_name', line.split('\t')[1])

#add GO terms
file = open('./data/gene_association.sgd')
lines = file.readlines()
file.close()
file = open('./data/gotermids')
gotermids = simplejson.load(file)
file.close()

for k in orfs.keys():
    setattr(orfs[k], 'goterms',{})

for line in lines:
    l = line.strip('\n').split('\t')
    if len(l) > 10:
        name = l[10].split('|')[0]
        if name in orfs.keys():
            goid = int(l[4].split(':')[1])
            orfs[name].goterms[l[8]] = gotermids[str(goid)]['name']

# add TFs and histone modifiers
file = open('./data/TFsJSON')
TFs = simplejson.load(file)
file.close()
histonefile = open('./data/histonemodifiersJSON')
histonemods = simplejson.load(histonefile)
histonefile.close()
for k in orfs.keys():
    if k in TFs:
        setattr(orfs[k], 'TF', True)
    else:
        setattr(orfs[k], 'TF', False)
    if hasattr(orfs[k], 'gene_name'):
        if orfs[k].gene_name in histonemods:
            setattr(orfs[k], 'histone_mod', True)
        else:
            setattr(orfs[k], 'histone_mod', False)
