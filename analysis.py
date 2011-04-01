import simplejson
import Paralogs
import ORFs
import globals

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
    Pars[(names[0], names[1])].alignment[names[0]] = dicpars[k][0]
    Pars[(names[0], names[1])].alignment[names[1]] = dicpars[k][1]

for p in plist:
    Paralogs.setalignindex(Pars, orfs, p)







