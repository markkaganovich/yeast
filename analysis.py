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

# add SYT positions 
for p in Pars.keys():
    if hasattr(Pars[p], 'alignindex'):
        setattr(Pars[p], 'sytpos', {})
        for o in Pars[p].orfs:
            if o in orfs.keys():
                Pars[p].sytpos[o] = (map(lambda x: Pars[p].alignindex[o][x], 
                            orfs[o].sytpos))

#within paralog seq conservation
# (what if there are more than two? - everything should be treated as a pair
# for neofunctionalization questions) -- right now everything is a pair
# do which amino acids are conserved, then how many
for p in Pars.keys():
    if hasattr(Pars[p], 'alignment'):
        o1 = Pars[p].orfs[0]
        o2 = Pars[p].orfs[1]
        seqcons = map(lambda x,y: x == y, Pars[p].alignment[o1], Pars[p].alignment[o2])
        sytconswithinSYT = list(set(Pars[p].sytpos[o1]) & set(Pars[p].sytpos[o2]))
        sytcons1 = ([x for x in Pars[p].sytpos[o1] if Pars[p].alignment[o1][x] ==
                Pars[p].alignment[o2][x]])
        sytcons2 = ([x for x in Pars[p].sytpos[o2] if Pars[p].alignment[o1][x] ==
                Pars[p].alignment[o2][x]])
        setattr(Pars[p], 'sytconswithinsyt', sytconswithinSYT)
        setattr(Pars[p], 'seqcons', seqcons)
        setattr(Pars[p], 'sytcons', sytcons1 + sytcons2)
        setattr(Pars[p], 'numseqcons', len(seqcons))
        setattr(Pars[p], 'numsytcons', len(Pars[p].sytcons))
        setattr(Pars[p], 'numsytconswithinsyt', len(sytconswithinSYT))


                








