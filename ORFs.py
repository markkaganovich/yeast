import simplejson
import globals

filehash = globals.json('./dbase/filehash')

class Orf:
    def __init__(self, name):
        self.name = name
	
    def getSYTpos(self):
        if hasattr(self, 'seq'):
	        sytpos = ([i for i, x in enumerate(self.seq) if x == 'S' or x == 's' or x == 'T' or x == 't' or x == 'Y' or x == 'y'])
	        setattr(self, 'sytpos', sytpos)

def addtoall(objects, geneset, attrfun, attrhash):
    for g in geneset:
        attrfun(objects[g], g, attrhash)

def getWparalogs(orf, wapinskiparalogskey):
    paralogs = []
    events = []
    name = orf.name
    wapinskiparalogs = globals.json(filehash[wapinskiparalogskey], globals.datasource)
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

def getphosphosites(orf, phosphositesfilekey):
    name = orf.name
    data = globals.json(filehash[phosphositesfilekey], globals.datasource)
    if name in data.keys():
        properlyindexed = map(lambda x: x-1, data[name])
        setattr(orf, 'phosphosite', properlyindexed)

def getalign(orf, filekey):
    name = orf.name
    data = globals.json(filehash[filekey], globals.datasource)
    if name in data.keys():
        setattr(orf, 'kwalalign', data[name])

def getotheralign(orf, *args):
    name = orf.name
    species = args[0]
    data = globals.json(filehash[species], globals.datasource)
    if name in data.keys():
        setattr(orf, species, data[name])

def getmultiplealign(orf, filekey):
    name = orf.name
    data = globals.json(filehash[filekey], globals.datasource)
    if name in data.keys():
        setattr(orf, 'multalign', data[name])

def getorthos(orf, filekey):
    orthohash = {}
    name = orf.name
    allorthologs = globals.json(filehash[filekey], globals.datasource)
    for key in allorthologs.keys():
        if name in allorthologs[key].keys():
            orthohash[key] =allorthologs[key][name]       
            setattr(orf, 'orthos', orthohash)
           
def getseq(orf, filekey):
    data = globals.json(filehash[filekey], globals.datasource)
    if orf.name in data.keys():
	    setattr(orf, 'seq', data[orf.name]) 


