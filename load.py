import ORFs
import Paralogs
import globals
import A
import iA
import iParA

def init(objects, clas, mod, types):
    g = getattr(A, types)
    data = globals.json(getattr(g, 'data'), globals.datasource)
    for t in data:
        objects[t] = clas(t)

def attrs(objects, internals, mod,  *attr):
    for a in attr:
        print a
        g = getattr(A, a)
        if hasattr(g, 'data'):
            data = globals.json(getattr(g, 'data'), globals.datasource)
        for o in objects.keys():
            if not hasattr(g, 'data'):
                if hasattr(g, 'args'):
                    args = getattr(g,'args')
                    mod.__getattribute__(getattr(g, 'fun'))(objects[o], args)
                else:
                    mod.__getattribute__(getattr(g, 'fun'))(objects[o])
            else:
                if hasattr(g, 'args'):
                    args = getattr(g,'args')
                    mod.__getattribute__(getattr(g, 'fun'))(objects[o], data, args)
                else:
                    mod.__getattribute__(getattr(g, 'fun'))(objects[o], data)
            intermethods = [i for i in internals.__dict__.keys() if not i.startswith('_')]
            apply(internal, [objects[o], internals] + sorted(intermethods))

def internal(obj,  internals, *internalattr):
    for ia in internalattr:
        g = getattr(internals, ia)
        if hasattr(g, 'args'):
            obj.__class__.__dict__[getattr(g, 'fun')](obj, getattr(g, 'args'))
        else:
            obj.__class__.__dict__[getattr(g, 'fun')](obj)

if __name__ == '__main__':
    orfs = {}
    init(orfs, ORFs.Orf, ORFs, 'init')
    attrs(orfs, iA, ORFs, 'orthologs', 'seq', 'kwalalign', 'phosphosites')

    pars = {}
    Paralogs.initall(pars, orfs, 'plistW')
    attrs(pars, iParA, Paralogs, 'paralign')
    attrs(pars, iParA, Paralogs,  'alignmentstats')
    Paralogs.calcdivergence(pars, 'kwalalignseqposcons', 'kwalalignsytposcons',
'kwalalignphosphositescons') 


