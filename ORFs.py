import simplejson

class Orf:
    def __init__(self, name):
        self.name = name

def getWparalogs(orf, name, wapinskiparalogs):
    paralogs = {}
    for event in wapinskiparalogs.keys():
        for k,v in wapinskiparalogs[event].items():
            if k == name:
                paralogs[name] = v
            if v == name:
                paralogs[name] = k

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

