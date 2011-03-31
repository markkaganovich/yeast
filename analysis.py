import simplejson
import Paralogs
import ORFs
'''
get wapinski paralogs, convert them to tupules, and add them as
attributes to the Paralog objects
''' 
file = open('./data/wapinskiparalogsbyevent')
paralogs = simplejson.load(file)
del paralogs['all']
file.close()

paraloglist = []
for key in paralogs.keys():
    for k,v in paralogs[key].items():
        paraloglist.append((k,v))

Pars = {}
Paralogs.initall(Pars, paraloglist)


