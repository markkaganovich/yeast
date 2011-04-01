import simplejson

class Orf:
    def __init__(self, name):
        self.name = name

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
    setattr(orf, 'phosphosite', phosphosites[name])

def getalign(orf, name, hash):
    if name in hash.keys():
        setattr(orf, 'kwalalign', hash[name])

def getmultiplealign(orf, name, hash):
    if name in hash.keys():
        setattr(orf, 'multalign', hash[name])

def getorthos(orf, name, allorthologs):
    orthohash = {}
    for key in allorthologs.keys():
        if name in allorthologs[key].keys():
            orthohash[key] =allorthologs[key][name]       
            setattr(orf, 'orthos', orthohash)

if __name__ == "__main__":
    file = open('./data/scergenes')
    geneset = simplejson.load(file)
    file.close()

    Orfs= {}
    file = open('./data/ScerKwalAlignment')
    scerkwalalign = simplejson.load(file)
    file.close()
    
    file = open('./data/speciesalign')
    speciesalign = simplejson.load(file)
    file.close()
     
    file = open('./data/KellisOrthos')
    orthologs = simplejson.load(file)
    file.close()
    
    file = open('./data/orthophosphosites')
    phosphosites = simplejson.load(file)
    file.close()
    
    file = open('./data/wapinskiparalogsbyevent')
    paralogs = simplejson.load(file)
    del paralogs['all']
    file.close()

def initall(objects, geneset):
    for g in geneset:
        objects[g] = Orf(g)

def addtoall(objects, geneset, attrfun, attrhash):
    for g in geneset:
        attrfun(objects[g], g, attrhash)

'''
        getalign(Orfs[g], g, scerkwalalign)
        getmultiplealign(Orfs[g], g, speciesalign)
        getorthos(Orfs[g], g, orthologs)
        getphosphosites(Orfs[g], g, phosphosites)
'''

'''
file = open('./data/sequences/Scer.fasta')
lines = file.readlines()
file.close()

for line in lines:
    if line.startswith('>'):
        o = orfs[line.strip('\n').split('>')[1]]
    else:
        setattr(o, 'seq', line.strip('\n'))
'''
