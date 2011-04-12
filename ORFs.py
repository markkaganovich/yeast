import simplejson

class Orf:
    def __init__(self, name):
        self.name = name

def initall(objects, geneset):
    for g in geneset:
        objects[g] = Orf(g)

def addtoall(objects, geneset, attrfun, attrhash):
    for g in geneset:
        attrfun(objects[g], g, attrhash)

def getWparalogs(orf, name, wapinskiparalogs):
    paralogs = []
    events = []
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

def getphosphosites(orf, name, phosphosites):
    if name in phosphosites.keys():
        properlyindexed = map(lambda x: x-1, phosphosites[name])
        setattr(orf, 'phosphosite', properlyindexed)

def getalign(orf, name, hash):
    if name in hash.keys():
        setattr(orf, 'kwalalign', hash[name])

def getotheralign(orf, name, hash, speciesalign):
    if name in hash.keys():
        setattr(orf, speciesalign, hash[name])


def getmultiplealign(orf, name, hash):
    if name in hash.keys():
        setattr(orf, 'multalign', hash[name])

def getorthos(orf, name, allorthologs):
    orthohash = {}
    for key in allorthologs.keys():
        if name in allorthologs[key].keys():
            orthohash[key] =allorthologs[key][name]       
            setattr(orf, 'orthos', orthohash)
           
def getSYTpos(orf):
    if hasattr(orf, 'seq'):
        sytpos = ([i for i, x in enumerate(orf.seq) if x == 'S' or x == 's' or 
                x== 'T' or x == 't' or x == 'Y' or x == 'y'])
        setattr(orf, 'sytpos', sytpos)


if __name__ == "__main__":
    file = open('./data/scergenes')
    geneset = simplejson.load(file)
    file.close()
    file = open('./data/ScerKwalAlignment')
    scerkwalalign = simplejson.load(file)
    file.close()
    file = open('./data/speciesalign')
    speciesalign = simplejson.load(file)
    file.close()
    file = open('./data/KellisOrthos')
    orthologs = simplejson.load(file)
    file.close()
    file = open('./data/gersteinphosphositefile')
    phosphosites = simplejson.load(file)
    file.close()
    file = open('./data/Scasalign')
    scasalign = simplejson.load(file)
    file.close()
    file = open('./data/Calbalign')
    calbalign = simplejson.load(file)
    file.close()
    file = open('./data/Spomalign')
    spomalign = simplejson.load(file)
    file.close()
    '''
    file = open('./data/orthophosphosites')
    phosphosites = simplejson.load(file)
    file.close()
    '''
    file = open('./paralogslist')
    paralogslist = simplejson.load(file)
    file.close()

    orfs={}
    initall(orfs, geneset)
    addtoall(orfs, geneset, getalign, scerkwalalign)
    addtoall(orfs, geneset, getmultiplealign, speciesalign)
    addtoall(orfs, geneset, getorthos, orthologs)
    addtoall(orfs, geneset, getphosphosites, phosphosites)
    addtoall(orfs, geneset, getWparalogs, paralogslist)
    
    map(lambda x: getSYTpos(orfs[x]), geneset)
    map(lambda x: getotheralign(orfs[x], x, scasalign, 'scasalign'), geneset)
    map(lambda x: getotheralign(orfs[x], x, calbalign, 'calbalign'), geneset)
    map(lambda x: getotheralign(orfs[x], x, spomalign, 'spomalign'), geneset)
    
    file = open('./data/sequences/Scer.fasta')
    lines = file.readlines()
    file.close()

    for line in lines:
        if line.startswith('>'):
            o = orfs[line.strip('\n').split('>')[1]]
        else:
            setattr(o, 'seq', line.strip('\n'))

