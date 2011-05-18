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

def turnover(alignment, index, ref, extent, binsize, aaset):
    turnovers = 0
    extentalign = alignment[extent]
    refalign = alignment[ref]
    topp = 0
    for i in aaset:
        ai = index[i]
        left = leftbinsalign(ai,  extentalign)
        right = rightbinsalign(ai,  extentalign)
        if refalign[ai] != extentalign[ai]:
            topp = topp+1
            indsl = left[0:binsize/2]
            indsr = right[0:binsize/2]
            tl = countturnovers(indsl, refalign, extentalign, ai)
            tr = countturnovers(indsr, refalign, extentalign, ai)
            if tl+tr >= 1:
                turnovers = turnovers + 1
    turndic = {}
    turndic['turnovers'] = turnovers
    turndic['turnopps'] = topp

def countturnovers(inds, refalign, extentalign, pos):
    for i in inds:
        if extentalign[i] == refalign[pos]:
            if refalign[i] != refalign[pos]:
                return 1
    return 0

def leftbinsalign(pos, alignment):
    left = [i for i,x in enumerate(alignment[0:pos]) if x != '-']
    left.reverse()
    return left

def rightbinsalign(pos, alignment):
    right = [i+pos+1 for i,x in enumerate(alignment[pos+1:]) if x != '-']
    return right
 


'''
clean up Pars -- get rid of orfs  not in orfs.keys(), move this to 
plist creation code from last week
'''
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
'''




                



