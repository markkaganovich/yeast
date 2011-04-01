import simplejson
import Paralogs

'''
get wapinski paralogs, convert them to tupules, and add them as
attributes to the Paralog objects
''' 

file = open('./data/paralogswapinskilistwithevents')
paraloglines = simplejson.load(file)
file.close()

paraloglist = map(lambda x: x.strip('\n').split(' '), paraloglines)
plist = [(x[0], x[1]) for x in paraloglist if len(x) > 2] 
'''
file = open('./paralogswapinskilist')
plist = simplejson.load(file)
file.close()
'''

Pars = {}
Paralogs.initall(Pars, plist)

file = open('./data/DICparsalignment')    
dicpars  = simplejson.load(file)
file.close()

for k in dicpars.keys():
    names = k.split('_')
    Pars[(names[0], names[1])].alignment = {}
    Pars[(names[0], names[1])].alignment[names[0]] = dicpars[k][0]
    Pars[(names[0], names[1])].alignment[names[1]] = dicpars[k][1]

for p in plist:
    Paralogs.setalignindex(Pars, orfs, p)

for p in Pars.keys():
    if hasattr(Pars[p], 'alignindex'):
        setattr(Pars[p], 'phosphosites', {})
        for o in Pars[p].orfs:
            if o in orfs.keys() and hasattr(orfs[o], 'phosphosite'):
                f = (filter(lambda x: x in Pars[p].alignindex[o].keys(), 
                            orfs[o].phosphosite))
                Pars[p].phosphosites[o] = (map(lambda x: Pars[p].alignindex[o][x],
                            f))






