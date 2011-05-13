import ORFs
import globals

def init(objects, clas, types):
    data = globals.json(filehash[types], globals.datasource)
    for t in data:
        objects[t] = clas(t)

def attrs(objects, *attr):
    for a in attr:
        data = globals.json(filehash[attrmethods[a]['data']], globals.datasource)
        for o in objects.keys():
            if 'args' in attrmethods[a].keys():
                args = attrmethods[a]['args']
                ORFs.__getattribute__(attrmethods[a]['fun'])(objects[o], data, args)
            else:
                ORFs.__getattribute__(attrmethods[a]['fun'])(objects[o], data)
            apply(internal, [objects[o]] + sorted(intermethods.keys()))

def internal(obj, *internalattr):
    for ia in internalattr:
        if 'args' in intermethods[ia].keys():
            obj.__class__.__dict__[intermethods[ia]['fun']](obj, intermethods[ia]['args'])
        else:
            obj.__class__.__dict__[intermethods[ia]['fun']](obj)
            
filehash = globals.json('./dbase/filehash')
attrmethods = globals.json('./dbase/attrmethods')
intermethods = globals.json('./dbase/intermethods')

if __name__ == '__main__':
    orfs = {}
    init(orfs, ORFs.Orf, 'geneset')




