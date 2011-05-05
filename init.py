import ORFs
import globals


def init(objects, clas, types):
	data  = globals.json(filehash[types], globals.datasource)
	for t in data:
		objects[t] = clas(t)
		
def attrs(objects, *attr):
	for a in attr:
		data = globals.json(filehash[a], globals.datasource)
		for o in objects.keys():
			ORFs.__getattribute__(attrmethods[a])(objects[o], o, data)

def internal(objects, *internalattr):
	for ia in internalattr:
		for o in objects.keys():
			objects[o].__class__.__dict__[intermethods[ia]]()

filehash = globals.json('./dbase/filehash')	
attrmethods = globals.json('./dbase/attrmethods')   
intermethods = globals.json('./dbase/intermethods') 
