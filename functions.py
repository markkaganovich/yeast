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

def getoneisTF(pars, attr):
    oneisTF = filter(lambda x: (hasattr(pars[x].orfs[0], attr) or 
hasattr(pars[x].orfs[1], attr)) and (getattr(pars[x].orfs[0], attr) == True or
getattr(pars[x].orfs[1], attr) == True), pars)
    return oneisTF

def getoneisWGDTF(pars, TFtype):
    oneisTF = getoneisTF(pars, TFtype)
    oneisWGDTF = filter(lambda x: pars[x].Wevent == 'WGD', oneisTF)
    return oneisWGDTF

def getoneisSSDTF(pars, TFtype):
    oneisTF = getoneisTF(pars, TFtype)
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

def diffgroups(pars, orfs, tfpairs, listofgroups, psitetype):
    groupmems = []
    tfmems = []
    for p in tfpairs:
        if hasattr(orfs[p[0]], 'genename') and hasattr(orfs[p[1]], 'genename'):
            a = map(lambda x: orfs[p[0]].genename in x, listofgroups)
            b = map(lambda x: orfs[p[1]].genename in x, listofgroups)
            print a
            print b
            gm = []
            if True in a and True in b:
                gm.append(a.index(True))
                gm.append(b.index(True))
                groupmems.append(gm)    
                tfmems.append([p[0], p[1]]) 
    diff = []
    difftfs = []
    for i,g in enumerate(groupmems):
        if g[0] - g[1] != 0:
            diff.append([g[0],g[1]]) 
            difftfs.append(tfmems[i])        
    
    difftfspsitessum = []
    difftfspsitesdiff = []
    alltfspsitessum = []
    alltfspsitesdiff = []
    for df in difftfs:
        difftfspsitessum.append(sum(map(lambda x: len(x),
getattr(pars[(df[0],df[1])], psitetype).values()))) 

    for a in tfpairs:
        alltfspsitessum.append(sum(map(lambda x: len(x), getattr(pars[(a[0],
a[1])], psitetype).values()))) 
           

    return [difftfspsitessum, alltfspsitessum]

def parselandry(lines, disorderedpsites):
    for l in lines:
        line = l.split('\t')
        if int(line[6]) == 1:# and (int(line[8]) !=1 and int(line[9]) != 1):
            if str(line[1]) in disorderedpsites.keys():
                disorderedpsites[str(line[1])].append(int(line[2]))
        
    return disorderedpsites


def calcpsiteoverlap(gersteinpsites, allpsites):
    total = 0
    for o in gersteinpsites.keys():
        if o in allpsites.keys():
            g = gersteinpsites[o]
            a = allpsites[o]
            total = total+len(filter(lambda x: x in g, a))

    return total
def allpsitedivs(pars):
    psitedivs = []
    for par in pars.keys():
        if hasattr(pars[par], 'psitediv'):
            psitedivs.append(pars[par].psitediv)
    return psitedivs

def phosphosdiv(pars, sets, psitetype):
    pdivs = []
    inds = []
    totalnum = []
    totals = []
    notconserveds = []
    ntotals = []
    conserved =[]
    constotals =[]
    numseqcons = []
    numseqconsrate= []
    totalnumdiff = []
    pdivs2 = []
    seqdivswpsites = []
    seqdivsnopsites = []
    eventsp = []
    eventsnop = []
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
            pdivs2.append(pdiv1+pdiv2)
            seqdiv1 = o1.kwalalignseqposcons
            seqdiv2 = o2.kwalalignseqposcons
            seqdivswpsites.append(abs(seqdiv1 - seqdiv2))
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
            eventsp.append(pars[s].Wevent)
        else:
            eventsnop.append(pars[s].Wevent)
            if hasattr(o1, 'kwalalignseqposcons') and hasattr(o2, 'kwalalignseqposcons'):
                seqdivsnopsites.append(abs(o1.kwalalignseqposcons - o2.kwalalignseqposcons)) 
        if hasattr(o1, psitetype):
            p1 = len(getattr(o1, psitetype))
        else:
            p1 = 0
        if hasattr(o2, psitetype):
            p2 = len(getattr(o2, psitetype))
        else:
            p2 = 0
        totalnum.append(p1+p2)
        totalnumdiff.append(abs(p1-p2))
        totals.append(totaldiff)
        ntotals.append(ntotaldiffs)
        constotals.append(constotaldiffs)
        conserved.append(cons)
        notconserveds.append(notcons)
        numseqcons.append(pars[s].numseqcons)
        length = len(o1.seq)+len(o2.seq)
        numseqconsrate.append(float(pars[s].numseqcons)/(length/2)) 
    return [eventsp, eventsnop, pdivs, pdivs2, seqdivsnopsites, seqdivswpsites, inds, totalnum, totals, notconserveds, ntotals, conserved,
constotals, numseqcons, numseqconsrate, totalnumdiff]


def calcevents(eventslist, pars, orfs):
    
    genes = [[]] * len(eventslist)
    pgenes = [[]] * len(eventslist)
    TFgenes = [[]] * len(eventslist)
    bulykTFs = [[]] * len(eventslist)
    pTFs = [[]] * len(eventslist)    
    TFpsites = [[]] * len(eventslist)
    bulykTFpsites = [[]] * len(eventslist)
    TFpgenes = [[]] * len(eventslist)
    bulykpgenes = [[]] * len(eventslist)
    psites = [[]] * len(eventslist)

    for i,e in enumerate(eventslist):
        genes[i] = filter(lambda x: pars[x].Wevent == e, pars)
        geneslist = []
        map(lambda x: geneslist.extend([x[0], x[1]]), genes[i])
        genes[i] = geneslist
        pgenes[i] = filter(lambda x: hasattr(orfs[x], 'phosphosite'), geneslist)
        psites[i] = map(lambda x: len(orfs[x].phosphosite), pgenes[i])
        TFgenes[i] = filter(lambda x: orfs[x].TF == True, geneslist)   
        TFpgenes[i] = filter(lambda x: hasattr(orfs[x], 'phosphosite'), TFgenes[i])
        TFpsites[i] = map(lambda x: len(orfs[x].phosphosite), TFpgenes[i])
        bulykTFs[i] = filter(lambda x: orfs[x].bulykTF == True, geneslist)
        bulykpgenes[i] = filter(lambda x: hasattr(orfs[x], 'phosphosite'), bulykTFs[i])  
        bulykTFpsites[i] = map(lambda x: len(orfs[x].phosphosite), bulykpgenes[i])
   
    return [genes, pgenes, psites, TFgenes, TFpgenes, TFpsites, bulykTFs, bulykpgenes,
bulykTFpsites]    

def getnondupl(orfs, pars):

    duplgenes = []
    map(lambda x: duplgenes.extend([x[0], x[1]]), pars)
    nondupl = filter(lambda x: x not in duplgenes, orfs.keys())

    psites = []
    TFs = []
    TFpsites = []

    pgenes = filter(lambda x: hasattr(orfs[x], 'phosphosite'), nondupl)
    psites = map(lambda x: len(orfs[x].phosphosite), pgenes)
    TFs = filter(lambda x: orfs[x].TF == True, nondupl)
    pTFs = filter(lambda x: hasattr(orfs[x], 'phosphosite'), TFs)
    TFpsites = map(lambda x: len(orfs[x].phosphosite), pTFs)    


    return [duplgenes, nondupl, psites, TFs, TFpsites]
def calcTFgroupspsites(group, orfs, psitetype):
    grouppsites = []
    for i in group:
        for o in orfs:
                if hasattr(orfs[o], 'genename') and orfs[o].genename == i:
                    if hasattr(orfs[o], psitetype):
                        grouppsites.append(len(getattr(orfs[o], psitetype)))
                    else:
                        grouppsites.append(0)
    return grouppsites

def orfpsites(orfs, group):
    psitenum = []
    psitecons = []
    genesused = []
    l = 0
    for i in group:
        for o in orfs:
            if hasattr(orfs[o], 'genename') and orfs[o].genename == i:        
                l = l+1       
                if hasattr(orfs[o], 'phosphosite') and hasattr(orfs[o],
'kwalalignphosphositecons'):
                    psitenum.append(len(orfs[o].phosphosite))
                    psitecons.append(orfs[o].kwalalignphosphositecons)
                    genesused.append(i)
                else:
                    psitenum.append(0)
    print l                
    return [psitenum, psitecons, genesused]

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


    

        
