import simplejson
import globals

filehash = globals.json('Paralogsfilehash', './dbase/')

class paralogs:
    def __init__(self, orfs):
        self.orfs = orfs

    def getalignindex(self):
        if hasattr(self, 'alignment') and getattr(self, 'alignment').keys().__len__() == 2:
            print 'alignindex:' + str(self.orfs[0].orfname) + str(self.orfs[1].orfname) 
            alignindex = {}
            dashes = 0
            i = 0
            alignseq = getattr(self, 'alignment')
            alignindex[self.orfs[0].orfname] = {}
            while i < len(self.orfs[0].seq):
                if alignseq[self.orfs[0].orfname][i+dashes] == '-':
                    dashes = dashes + 1
                else:
                    alignindex[self.orfs[0].orfname][i] = dashes+i
                    i = i+1
            alignindex[self.orfs[1].orfname] = {}
            i = 0
            dashes = 0
            while i < len(self.orfs[1].seq):
                if alignseq[self.orfs[1].orfname][i+dashes] == '-':
                    dashes = dashes + 1
                else:
                    alignindex[self.orfs[1].orfname][i] = dashes+i
                    i = i+1
            setattr(self, 'alignindex', alignindex)
    
    def getphosphosites(self):
        setattr(self, 'phosphosites', {})
        print 'phosphosite:'+ str(self.orfs[0].orfname)
        print 'phosphosite:' + self.orfs[1].orfname
        for i in range(0,len(self.orfs)):
            if hasattr(self.orfs[i], 'phosphosite') and hasattr(self,
'alignindex'):
                self.phosphosites[self.orfs[i].orfname] = filter(lambda x: self.alignindex[self.orfs[i].orfname][x], self.orfs[i].phosphosite)

    def getsyt(self):
        if hasattr(self, 'alignindex'):
            sytpos = {}
            setattr(self, 'sytpos', sytpos)
            for i in self.orfs:
                self.sytpos[i.orfname] = (map(lambda x: self.alignindex[i.orfname][x], i.sytpos))

def alignmentstats(self):
    if hasattr(self, 'alignment'):
        o1 = self.orfs[0]
        o2 = self.orfs[1]
        seqcons = map(lambda x,y: x == y, self.alignment[o1.orfname], self.alignment[o2.orfname])
        sytconswithinSYT = list(set(self.sytpos[o1.orfname]) & set(self.sytpos[o2.orfname]))
        sytcons1 = ([x for x in self.sytpos[o1.orfname] if self.alignment[o1.orfname][x] == self.alignment[o2.orfname][x]])
        sytcons2 = ([x for x in self.sytpos[o2.orfname] if self.alignment[o1.orfname][x] == self.alignment[o2.orfname][x]])
        setattr(self, 'sytconswithinsyt', sytconswithinSYT)
        setattr(self, 'seqcons', seqcons)
        setattr(self, 'sytcons', sytcons1 + sytcons2)
        setattr(self, 'numseqcons', sum(seqcons))
        setattr(self, 'sumseqlength', len(o1.seq)+len(o2.seq))
        setattr(self, 'seqdiff', 1-(self.numseqcons/float(self.sumseqlength/2)))
        setattr(self, 'numsytcons', len(list(set(self.sytcons))))
        setattr(self, 'numsytconswithinsyt', len(list(set(sytconswithinSYT))))
        setattr(self, 'avgsyt', (len(self.sytpos[o1.orfname])+len(self.sytpos[o2.orfname])/2))
        setattr(self, 'sytdiff', 1-(float(self.numsytcons)/self.avgsyt))

def phosphostats(self):
    o1 = self.orfs[0]
    o2 = self.orfs[1]
    if hasattr(self, 'phosphosites'):
        if len(self.phosphosites.keys()) != 2:
            setattr(self, 'phosphositecons', 0)
        else:
            phosphositecons = len(list(set(self.phosphosites[o1.orfname]) &
set(self.phosphosites[o2.orfname])))
            (setattr(self, 'phosphositecons', phosphositecons/
                     (float((len(self.phosphosites[o1.orfname]) +
len(self.phosphosites[o2.orfname])))/2)))

def getevents(self, data):
    key = self.orfs[0].orfname + '_' + self.orfs[1].orfname
    setattr(self, 'Wevent', data[key])

def initall(objects, orfs, parset):
    plist = globals.json(filehash[parset])
    for par in plist:
        objects[(par[0], par[1])] = paralogs((orfs[par[0]], orfs[par[1]]))

def calcdivergence(objects, seqcons, sytcons, psitecons):
    for p in objects.keys():
        par = objects[p]
        o1 = par.orfs[0]
        o2 = par.orfs[1]
        if o1.orthos['Kwal'] == 'NONE' or o2.orthos['Kwal'] == 'NONE':
            continue
        setattr(par, 'minorfseqdiv', min(1-getattr(o1,seqcons), 1-getattr(o2,seqcons)))
        setattr(par, 'maxorfseqdiv', max(1-getattr(o1,seqcons), 1-getattr(o2,seqcons)))
        setattr(par, 'minorfsytdiv', min(1-getattr(o1,sytcons), 1-getattr(o2,sytcons)))
        setattr(par, 'maxorfsytdiv', max(1-getattr(o1,sytcons), 1-getattr(o2,sytcons)))
        if hasattr(o1, psitecons) and hasattr(o2, psitecons):
            setattr(par, 'minorfphosphodiv', min(1-getattr(o1,psitecons), 1-getattr(o2,psitecons)))
            setattr(par, 'maxorfphosphodiv', max(1-getattr(o1,psitecons),1-getattr(o2,psitecons)))


def getalign(par, alignments):
    alignment = {}
    name1 = par.orfs[0].orfname + '_' + par.orfs[1].orfname
    name2 = par.orfs[1].orfname + '_' + par.orfs[0].orfname
    if name1 in alignments.keys():
        name = name1
    elif name2 in alignments.keys():
        name = name2
    else:
        print par.orfs[0].orfname
        print par.orfs[1].orfname
        return
    alignment = alignments[name]
    setattr(par, 'alignment', alignment)

def reminstances(objects, attr):
    for p in objects.keys():
        if not hasattr(objects[p], attr):
            del objects[p]




        
        


