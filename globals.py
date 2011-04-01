import simplejson

def addtoall(objects, geneset, attrfun, attrhash):
    for g in geneset:
        attrfun(objects[g], g, attrhash)

def writebyparts(serializeddict, filepath):
    for i in range(0,10):
        newdic = {}
        start = i*(len(serializeddict) / 10)
        file = open(filepath+str(i),'w')
        for j in range(0, len(serializeddict)/10):
            key = serializeddict.keys()[start+j]
            newdic[key] = serializeddict[key]
        simplejson.dump(newdic, file)

def readbyparts(filepath):
    dic = {}
    for i in range(0,10):
        file = open(filepath+str(i))
        dic.update(simplejson.load(file))
        file.close()

    return dic




def serialize(dictobj):
    newdict = {}
    for key in dictobj.keys():
        newdict[key] = dictobj[key].__dict__

    return newdict



