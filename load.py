import ORFs
import Paralogs
import globals
import A
import iA

def init(objects, clas, mod, types):
    g = getattr(A, types)
    data = globals.json(getattr(g, 'data'), globals.datasource)
    for t in data:
        objects[t] = clas(t)

def attrs(objects, mod,  *attr):
    for a in attr:
        print a
        g = getattr(A, a)
        data = globals.json(getattr(g, 'data'), globals.datasource)
        for o in objects.keys():
            if hasattr(g, 'args'):
                args = getattr(g,'args')
                ORFs.__getattribute__(getattr(g, 'fun'))(objects[o], data, args)
            else:
                ORFs.__getattribute__(getattr(g, 'fun'))(objects[o], data)
            intermethods = [i for i in iA.__dict__.keys() if not
i.startswith('_')]
            apply(internal, [objects[o]] + sorted(intermethods))

def internal(obj,  *internalattr):
    for ia in internalattr:
        g = getattr(iA, ia)
        print g
        if hasattr(g, 'args'):
            obj.__class__.__dict__[getattr(g, 'fun')](obj, getattr(g, 'args'))
        else:
            obj.__class__.__dict__[getattr(g, 'fun')](obj)

if __name__ == '__main__':
    orfs = {}
    init(orfs, ORFs.Orf, ORFs, 'init')
    attrs(orfs, ORFs, 'seq')

    pars = {}
    Paralogs.initall(pars, orfs, 'plistW')

    


