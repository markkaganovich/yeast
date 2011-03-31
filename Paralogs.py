import simplejson

#parset is list of tupules
def initall(objects, parset):
   for par in parset:
       objects[par] = paralogs(par)

class paralogs:
    def __init__(self, *orfs):
        self.orfs = orfs


        
        


