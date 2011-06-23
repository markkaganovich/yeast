import globals
import numpy as np
'''
input bi-alignment, output indeces of alignment to find SVs
'''
def SVs(alignment, slicesize):
    SVs = {}
    SVspos = {}
    count = 0
    positions = []
    for s in alignment.keys():
        startstop = []
        for i,x in enumerate(alignment[s]):
            if x == '-':
                count = count + 1
            else:
                if count >= slicesize:
                    startstop.append((i - count, i))
                    positions.extend(range(i-count, i))
                count = 0
        SVspos[s] = startstop
        SVs = positions
    return SVs

def translatetoseqpos(object):
    if hasattr(object, 'kwalalignSVs'):
        seqinds = map(lambda x: reverse_lookup(object.kwalalignindex, x),
object.kwalalignSVs)
        setattr(object, 'SVseqinds', seqinds)

def addtranslateseqpos(objects):
    for o in objects:
        translatetoseqpos(objects[o])

def addtoobjects(objects, alignmentname, slicesize):
    for k in objects:
        o = objects[k]
        if hasattr(o, alignmentname):
            alignment = getattr(o, alignmentname) 
            svs = SVs(alignment, slicesize)
            setattr(o, str(alignmentname)+'SVs', svs)
 

# do total SV length between paralogs
# SV length correlates to TF divergence more than seq random sequence
# divergence?

# SV positions and SV + phosphosite positions intersected

# are there paralogs that both have overlapping SVs wrt Orthologs in KWAL?

def SVintersect(self):
    o1 = self.orfs[0]
    o2 = self.orfs[1]
    if hasattr(o1, 'SVseqinds') and hasattr(o2, 'SVseqinds'):
        SVpos1 = getattr(o1, 'SVseqinds')
        SVpos2 = getattr(o2, 'SVseqinds')
    intx = list(set(SVpos1) and set(SVpos2))
    setattr(self, 'SVintersect', intx.__len__())
    return intx.__len__()

def getWGD(pars):
    WGDpars = filter(lambda x: pars[x].Wevent == 'WGD', pars)
    return WGDpars

def getSSD(pars):
    SSDpars = filter(lambda x: pars[x].Wevent != 'WGD', pars)
    return SSDpars

def getbothparTFs(pars):
    bothparTF = filter(lambda x: hasattr(pars[x].orfs[0], 'TF') and
hasattr(pars[x].orfs[1], 'TF') and pars[x].orfs[0].TF == True and pars[x].orfs[1].TF == True, pars)
    return bothparTF

def getWGDparTFs(pars):
    bothparTFs = getbothparTFs(pars) 
    bothWGDparTFs = filter(lambda x: pars[x].Wevent == 'WGD', bothparTFs)
    return bothWGDparTFs

def getoneisTF(pars):
    oneisTF = filter(lambda x: (hasattr(pars[x].orfs[0], 'TF') or 
hasattr(pars[x].orfs[1], 'TF')) and (pars[x].orfs[0].TF == True or
pars[x].orfs[1].TF == True), pars)
    return oneisTF

def getoneisWGDTF(pars):
    oneisTF = getoneisTF(pars)
    oneisWGDTF = filter(lambda x: pars[x].Wevent == 'WGD', oneisTF)
    return oneisWGDTF

def getoneisSSDTF(pars):
    oneisTF = getoneisTF(pars)
    oneisSSDTF = filter(lambda x: pars[x].Wevent != 'WGD', oneisTF)
    return oneisSSDTF

def psiteVSexprd(pars,parTFs):
    psitedivs = []
    inds = []
    for i,par in enumerate(parTFs):
        if hasattr(pars[par], 'psitediv'): 
            psitedivs.append(pars[par].psitediv)
            inds.append(i)
            print pars[par].phosphosites
    return [psitedivs, inds]

def allpsitedivs(pars):
    psitedivs = []
    for par in pars.keys():
        if hasattr(pars[par], 'psitediv'):
            psitedivs.append(pars[par].psitediv)
    return psitedivs

def phosphosdiv(pars, sets):
    pdivs = []
    inds = []
    totalnum = []
    totals = []
    notconserveds = []
    ntotals = []
    conserved =[]
    constotals =[]
    for i,s in enumerate(sets):
        o1 = pars[s].orfs[0]
        o2 = pars[s].orfs[1]
        totaldiff = 0
        ntotaldiffs = 0
        constotaldiffs = 0
        cons = 0
        notcons = 0
        
        if hasattr(o1, 'kwalalignphosphositecons') and hasattr(o2, 'kwalalignphosphositecons'):
            pdiv1 = o1.kwalalignphosphositecons
            pdiv2 = o2.kwalalignphosphositecons
            pdivs.append(abs(pdiv1-pdiv2))
            inds.append(i)
            t = len(pars[s].phosphosites.values()[0]) + len(pars[s].phosphosites.values()[1])
            u = len(list(set(pars[s].phosphosites.values()[0]) & set(pars[s].phosphosites.values()[1])))
            totaldiff = t-u
            nc1 = map(lambda x: o1.kwalalignindex[x], o1.notconservedpsites)
            nc2 = map(lambda x: o2.kwalalignindex[x], o2.notconservedpsites)
            nt = len(nc1) + len(nc2)
            nu = len(list(set(nc1) & set(nc2)))
            notcons = len(o1.notconservedpsites)+len(o2.notconservedpsites)
            cons = len(o1.conservedpsites)+len(o2.conservedpsites)
            c1 = map(lambda x: o1.kwalalignindex[x], o1.conservedpsites)
            c2 = map(lambda x: o2.kwalalignindex[x], o2.conservedpsites)
            ct = len(c1) + len(c2)
            cu = len(list(set(c1) & set(c2)))
            constotaldiffs = ct - cu
            ntotaldiffs = nt-nu
        if hasattr(o1, 'phosphosite'):
            p1 = len(o1.phosphosite)
        else:
            p1 = 0
        if hasattr(o2, 'phosphosite'):
            p2 = len(o2.phosphosite)
        else:
            p2 = 0
        totalnum.append(p1+p2)
        totals.append(totaldiff)
        ntotals.append(ntotaldiffs)
        constotals.append(constotaldiffs)
        conserved.append(cons)
        notconserveds.append(notcons)
    return [pdivs, inds, totalnum, totals, notconserveds, ntotals, conserved,
constotals]
    
def orfphosphodiv(orfs, pair1, pair2):
    pdivs = []
    totalnum = []
    for i in range(0,len(pair1)):
        if pair1[i] not in orfs.keys() or pair2[i] not in orfs.keys():
            continue
        o1 = orfs[pair1[i]]
        o2 = orfs[pair2[i]]
        if hasattr(o1, 'kwalalignphosphositecons') and hasattr(o2, 'kwalalignphosphositecons'):
            pdiv1 = o1.kwalalignphosphositecons
            pdiv2 = o2.kwalalignphosphositecons
            pdivs.append(abs(pdiv1-pdiv2))
    if hasattr(o1, 'phosphosite'):
        p1 = len(o1.phosphosite)
    else:
        p1 = 0
    if hasattr(o2, 'phosphosite'):
        p2 = len(o2.phosphosite)
    else:
        p2 = 0
    totalnum.append(p1+ p2)
    return [pdivs, totalnum]


def simulation(orfs):
    genes = orfs.keys()
    pair = choose_without_replacement(len(genes), 2)
    [simpdiv, simtotalnum] = orfphosphodiv(orfs, [genes[pair[0]-1]],
[genes[pair[1]-1]])    
    return [genes[pair[0]-1], genes[pair[1]-1], simpdiv, simtotalnum]
    
def runsim(orfs):
    pair1 = []
    pair2 = []
    simpdivs = []
    simtotalnums = []
    for i in range(0,10000):
        a = simulation(orfs)
        simpdivs.append(a[2])
        simtotalnums.append(a[3])
        pair1.append(a[0])
        pair2.append(a[1])
    return [pair1, pair2, simpdivs, simtotalnums]


def choose_without_replacement(m,n,repeats=None):
   """Choose n nonnegative integers less than m without replacement

   Returns an array of shape n, or (n,repeats).
   """
   if repeats is None:
       r = 1
   else:
       r = repeats
   if n>m:
       raise ValueError, "Cannot find %d nonnegative integers less than %d" %  (n,m)
   if n>m/2:
       res = np.sort(np.random.rand(m,r).argsort(axis=0)[:n,:],axis=0)
   else:
       res = np.random.random_integers(m,size=(n,r))
       while True:
           res = np.sort(res,axis=0)
           w = np.nonzero(np.diff(res,axis=0)==0)
           nr = len(w[0])
           if nr==0:
               break
           res[w] = np.random.random_integers(m,size=nr)

   if repeats is None:
       return res[:,0]
   else:
       return res


def reverse_lookup(d, v):
    for k in d:
        if d[k] == v:
            return k


    

        
