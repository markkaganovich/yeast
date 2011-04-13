'''
check for turnover
    syt
    all amino acids
    p-sites

if the amino acid is mutated in the alignment, does it
appear within X amino acids (in the sequence being compared
to the REF).... and same thing but if that aa 'appears'.. i.e.
the ref is a different aa

do this for ortho-para comparison and para-para comparison
'''

syt = ['S', 'Y', 'T']
aaset = syt
binsize = 10
for p in Pars.keys():
    par = Pars[p]
    o1 = par.orfs[0]
    o2 = par.orfs[1]
    setattr(par, 'lost',{})
    for o in [orfs[o1],orfs[o2]]:
        par.lost[o] = 0
        if 'Scer' in o.kwalalign.keys():
            for i,x in enumerate(o.kwalalign['Scer']):
                if x in aaset and x!= o.kwalalign['Kwal'][i]:
                    slice = o.kwalalign['Kwal'].__getslice__(i-binsize/2,i+binsize/2)    
                    inds = [l for l in range(0, len(slice)) if slice[l] == x]
                    lost = filter(lambda y: o.kwalalign['Scer'][y] != x, inds)
                    if lost != []:
                        par.lost[o] = len(lost) + par.lost[o]
                    
'''
turnover in orfs
'''
syt = ['S', 'Y', 'T']
aaset = syt
binsize = 10

def calcturnoverorfs(orfs, binsize, aapos):
    for k in orfs.keys():
        o = orfs[k]
        aaset = getattr(o, aapos)
        if 'Scer' in o.kwalalign.keys():
            o.turnover[binsize] = 0
            for i,x in enumerate(o.kwalalign['Scer']):
                if i in aaset and x != o.kwalalign['Kwal'][i]:
                    slice = o.kwalalign['Kwal'].__getslice__(i-binsize/2,i+binsize/2)    
                    inds = [l for l in range(0, len(slice)) if slice[l] == x]
                    lost = filter(lambda y: o.kwalalign['Scer'][y] != x, inds)
                    if lost != []:
                        o.turnover[binsize] = len(lost) + o.turnover[binsize]

binsizes = range(0,300,10)
for b in binsizes:
    calcturnoverorfs(orfs, b, aaset)



'''
clean up Pars -- get rid of orfs  not in orfs.keys(), move this to 
plist creation code from last week
'''
allparsorfs = []
for p in Pars.keys():
    if p[0] not in orfs.keys():
        del Pars[p]
    elif p[1] not in orfs.keys():
        del Pars[p]
    else:
        allparsorfs.append(p[0])
        allparsorfs.append(p[1])



                



