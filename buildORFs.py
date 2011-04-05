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
