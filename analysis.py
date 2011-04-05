import simplejson
import Paralogs
'''
get wapinski paralogs, convert them to tupules, and add them as
attributes to the Paralog objects
''' 
file = open('./data/paralogswapinskilistfigS6')
paraloglines = simplejson.load(file)
file.close()

paraloglist = map(lambda x: x.strip('\n').split('\t'), paraloglines)
plist = [(x[0], x[1]) for x in paraloglist] 
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
        setattr(Pars[p], 'numseqcons', sum(seqcons))
        setattr(Pars[p], 'sumseqlength', len(orfs[o1].seq)+len(orfs[o2].seq))
        setattr(Pars[p], 'seqdiff', 1- (Pars[p].numseqcons/float(Pars[p].sumseqlength/2))) 

for p in Pars.keys():
    if hasattr(Pars[p], 'alignment'):
        o1 = Pars[p].orfs[0]
        o2 = Pars[p].orfs[1]
        setattr(Pars[p], 'numsytcons', len(list(set(Pars[p].sytcons))))
        setattr(Pars[p], 'numsytconswithinsyt', len(list(set(sytconswithinSYT))))
        setattr(Pars[p], 'avgsyt', (len(Pars[p].sytpos[o1])+len(Pars[p].sytpos[o2])/2))
        setattr(Pars[p], 'sytdiff', 1-(float(Pars[p].numsytcons)/Pars[p].avgsyt))

for p in Pars.keys():
    o1 = Pars[p].orfs[0]
    o2 = Pars[p].orfs[1]
    if hasattr(Pars[p], 'phosphosites'):
        if len(Pars[p].phosphosites.keys()) != 2:
            setattr(Pars[p], 'phosphositecons', 0)
        else:
            phosphositecons = len(list(set(Pars[p].phosphosites[o1]) & set(Pars[p].phosphosites[o2])))
            (setattr(Pars[p], 'phosphositecons', phosphositecons/
                     (float((len(Pars[p].phosphosites[o1]) + len(Pars[p].phosphosites[o2])))/2)))
        
        '''    
        if len(Pars[p].phosphosites.keys()) == 0:
            setattr(Pars[p], 'totalphosphoevents', 0)
        else:
            phosphos = 0
            for o in (o1,o2):
                if o in Pars[p].phosphosites.keys():
                    phosphos = phosphos + len(Pars[p].phosphosites[o])
        '''
# add events 
events = [x[2] for x in paraloglist if len(x) > 2] 
for i,p in enumerate(plist):
    setattr(Pars[p], 'event', events[i])
    
gotermdiv = [(x[3], x[4], x[5]) for x in paraloglist if len(x) > 5]
for i,p in enumerate(plist):
    setattr(Pars[p], 'go_component_div', gotermdiv[i][0])
    setattr(Pars[p], 'go_process_div', gotermdiv[i][1])
    setattr(Pars[p], 'go_function_div', gotermdiv[i][2])

txnmodulediv = [x[6] for x in paraloglist]
for i,p in enumerate(plist):
    setattr(Pars[p], 'txnmodule_div_wapin', txnmodulediv[i])





    










