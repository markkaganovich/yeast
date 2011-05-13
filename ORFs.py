import simplejson
import globals

filehash = globals.json('./dbase/filehash')

class Orf:
    def __init__(self, orfname):
        self.orfname = orfname
	
    def getSYTpos(self):
        if hasattr(self, 'seq'):
	        sytpos = ([i for i, x in enumerate(self.seq) if x == 'S' or x == 's' or 
            x == 'T' or x == 't' or x == 'Y' or x == 'y'])
	        setattr(self, 'sytpos', sytpos)

    def getalignindex(self, align):
        if hasattr(self, align) and getattr(self, align).keys().__len__() == 2:
            alignindex = {}
            dashes = 0
            i = 0
            alignseq = getattr(self, align)
            while i < len(self.seq):
                if alignseq['Scer'][i+dashes] == '-':
                    dashes = dashes + 1
                else:
                    alignindex[i] = dashes+i
                    i = i+1
            setattr(self, align+'index', alignindex)

    def seqconservation(self, args):
        align = args[0]
        aaposition = args[1]
        if (hasattr(self, align) and getattr(self, align).keys().__len__() == 2 
            and hasattr(self, aaposition)):
            alignseq = getattr(self, align)
            otherspecies = alignseq.keys().remove('Scer')
            index = map(lambda x: getattr(self, align+'index')[x], getattr(self, aaposition))
            cons = map(lambda x: alignseq[alignseq.keys()[0]][x] == alignseq[alignseq.keys()[1]][x], index)
            setattr(self, align+aaposition+'cons', sum(cons)/float(len(getattr(self, aaposition))))
        

def addtoall(objects, geneset, attrfun, attrhash):
    for g in geneset:
        attrfun(objects[g], g, attrhash)

def genename(orf, genenames):
    if orf.orfname in genenames.keys():
        orf.genename = genenames[orf.orfname]

def getWparalogs(orf, wapinskiparalogs):
    paralogs = []
    events = []
    name = orf.orfname
    for p in wapinskiparalogs:
        l = p.strip('\n').split(' ')
        if len(l) > 2:
            if l[0] == name:
                paralogs.append(l[1])
                events.append(l[2])
            if l[1] == name:
                paralogs.append(l[0])
                events.append(l[2])
    setattr(orf, 'paralogs', paralogs)  
    setattr(orf, 'duplevent', events)

def getphosphosites(orf, data):
    name = orf.orfname
    if name in data.keys():
        properlyindexed = map(lambda x: x-1, data[name])
        setattr(orf, 'phosphosite', properlyindexed)

def getalign(orf, data):
    name = orf.orfname
    if name in data.keys():
        setattr(orf, 'kwalalign', data[name])

def getotheralign(orf, data, *args):
    name = orf.orfname
    species = args[0]
    if name in data.keys():
        setattr(orf, species, data[name])

def getmultiplealign(orf, data):
    name = orf.orfname
    if name in data.keys():
        setattr(orf, 'multalign', data[name])

def getorthos(orf, allorthologs):
    orthohash = {}
    name = orf.orfname
    for key in allorthologs.keys():
        if name in allorthologs[key].keys():
            orthohash[key] =allorthologs[key][name]       
            setattr(orf, 'orthos', orthohash)
           
def getseq(orf, data):
    if orf.orfname in data.keys():
        setattr(orf, 'seq', data[orf.orfname])
        setattr(orf, 'seqpos', range(0,len(orf.seq)))

def hasproperty(orf, data, attr):
    if orf.orfname in data or (hasattr(orf, 'genename') and orf.genename in
data):
        setattr(orf, attr, True)
    else:
        setattr(orf, attr, False)

def getGOterms(orf, data):
    setattr(orf, 'goterms', data[orf.orfname])


    



