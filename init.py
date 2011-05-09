import ORFs
import globals


def init(objects, clas, types):
	data  = globals.json(filehash[types], globals.datasource)
	for t in data:
		objects[t] = clas(t)
		
def attrs(objects, *attr):
	for o in objects.keys():
		for a in attr:
			if 'args' in attrmethods[a].keys():
		       args = attrmethods[a]['args']
			   ORFs.__getattribute__(attrmethods[a]['fun'])(objects[o], args)
			apply(internal, [objects[o]] + intermethods.keys())

def internal(obj, *internalattr):
	for ia in internalattr:
			obj.__class__.__dict__[intermethods[ia]](obj)

filehash = globals.json('./dbase/filehash')	
attrmethods = globals.json('./dbase/attrmethods')   
intermethods = globals.json('./dbase/intermethods') 

if __name__ == '__main__':
	orfs ={}
	init(orfs, ORFs.Orf, 'geneset')
