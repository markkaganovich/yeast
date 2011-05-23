import simplejson
import globals

filehash = globals.json('Paralogsfilehash', './dbase/')
attrmethods = globals.json('Paralogsattrmethods', './dbase/')
intermethods = globals.json('Paralogsintermethods', './dbase/')

class paralogs:
    def __init__(self, orfs):
        self.orfs = orfs

    def getalignindex(self):
        if hasattr(self, 'alignment') and getattr(self, 'alignment').keys().__len__() == 2:
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


def initall(objects, orfs, parset):
    plist = globals.json(filehash[parset])
    for par in plist:
        objects[(par[0], par[1])] = paralogs((orfs[par[0]], orfs[par[1]]))

def getalign(par, alignments):
    alignment = {}
    alignment[par.orfs[0].orfname] = alignments[par.orfs[0].orfname+'_'+par.orfs[1].orfname][0]
    alignment[par.orfs[1].orfname] = alignments[par.orfs[0].orfname+'_'+par.orfs[1].orfname][1]
    setattr(par, 'alignment', alignment)




        
        


