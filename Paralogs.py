import simplejson

#parset is list of tupules
def initall(objects, parset):
   for par in parset:
       objects[par] = paralogs(par)

def setalignindex(objects, orfs, name):
    alignindex={}
    alignindex[name[0]]={}
    alignindex[name[1]] = {}
    if hasattr(objects[name], 'alignment'):
        for k in range(0,2):
            dashes = 0
            print name
            if name[k] in orfs.keys() and hasattr(orfs[name[k]], 'seq'):
                for i in range(0, len(orfs[name[k]].seq)):
                    if objects[name].alignment[name[k]][i+dashes] == '-':
                        dashes = dashes + 1
                    alignindex[name[k]][i] = dashes +i 
        setattr(objects[name], 'alignindex', alignindex)


class paralogs:
    def __init__(self, orfs):
        self.orfs = orfs



        
        


