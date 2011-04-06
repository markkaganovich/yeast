'''
yeast gene duplication results in new function, bitches

the existence of orfs and Pars is assumed
'''
import fisher
import numpy

'''
compare phosphosite divergence with sequence divergence
'''
divGOseqdiv = {}
divGOphosphodiv = {}
divGOsytdiv = {}
consGOseqdiv = {}
consGOphosphodiv = {}
consGOsytdiv = {}
types = ['P','C','F']
typesinPars = ['go_process_div','go_function_div','go_component_div']

wgdpar = filter(lambda x: Pars[x].event == 'WGD', Pars.keys())
ssd = [x for x in Pars.keys() if Pars[x].orfs[0] in orfs.keys() and Pars[x].orfs[1] in orfs.keys() and orfs[Pars[x].orfs[0]].orthos['Kwal'] == orfs[Pars[x].orfs[1]].orthos['Kwal'] and Pars[x].event != 'WGD']
ssdbefore = [x for x in Pars.keys() if Pars[x].orfs[0] in orfs.keys() and Pars[x].orfs[1] in orfs.keys() and orfs[Pars[x].orfs[0]].orthos['Calb'] == orfs[Pars[x].orfs[1]].orthos['Calb'] and Pars[x].event != 'WGD']

'''
ssd1 = [x for x in Pars.keys() if Pars[x].orfs[0] in orfs.keys() and Pars[x].orfs[1] in orfs.keys() and orfs[Pars[x].orfs[0]].orthos['Spar'] == orfs[Pars[x].orfs[1]].orthos['Spar'] and Pars[x].event != 'WGD']
ssd2 = [x for x in Pars.keys() if not x in ssd1 and Pars[x].orfs[0] in orfs.keys() and Pars[x].orfs[1] in orfs.keys() and orfs[Pars[x].orfs[0]].orthos['Sbay'] == orfs[Pars[x].orfs[1]].orthos['Sbay'] and Pars[x].event != 'WGD']
'''
ssd3 = [x for x in Pars.keys() if Pars[x].orfs[0] in orfs.keys() and Pars[x].orfs[1] in orfs.keys() and 'Scas' in orfs[Pars[x].orfs[0]].orthos.keys() and 'Scas' in orfs[Pars[x].orfs[1]].orthos.keys() and orfs[Pars[x].orfs[0]].orthos['Scas'] == orfs[Pars[x].orfs[1]].orthos['Scas'] and Pars[x].event != 'WGD']
ssd4 = [x for x in Pars.keys() if not x in ssd3 and Pars[x].orfs[0] in orfs.keys() and Pars[x].orfs[1] in orfs.keys() and orfs[Pars[x].orfs[0]].orthos['Kwal'] == orfs[Pars[x].orfs[1]].orthos['Kwal'] and Pars[x].event != 'WGD']
ssd5 = [x for x in Pars.keys() if not x in ssd4 and Pars[x].orfs[0] in orfs.keys() and Pars[x].orfs[1] in orfs.keys() and orfs[Pars[x].orfs[0]].orthos['Calb'] == orfs[Pars[x].orfs[1]].orthos['Calb'] and Pars[x].event != 'WGD']
ssd6 = [x for x in Pars.keys() if not x in ssd5 and Pars[x].orfs[0] in orfs.keys() and Pars[x].orfs[1] in orfs.keys() and orfs[Pars[x].orfs[0]].orthos['Spom'] == orfs[Pars[x].orfs[1]].orthos['Spom'] and Pars[x].event != 'WGD']

ssdunion = ssd3 +ssd4 + ssd5 + ssd6

type = 'go_process_div'
divGOmaxseqdiv = []
divGOmaxphosphodiv = [] 
divGOmaxsytdiv = [] 
divGOminseqdiv = [] 
divGOminphosphodiv = [] 
divGOminsytdiv = []
consGOmaxseqdiv = [] 
consGOmaxphosphodiv = []
consGOmaxsytdiv = [] 
consGOminseqdiv = [] 
consGOminphosphodiv = [] 
consGOminsytdiv = []
here =[]
for p in geneset:
    par = Pars[p]
    if hasattr(par, 'maxorfseqdiv'):
        here.append(p)
        if getattr(par, type) == str(1):
            divGOmaxseqdiv.append(par.maxorfseqdiv)
            divGOmaxsytdiv.append(par.maxorfsytdiv)
            if hasattr(par, 'maxorfphosphodiv'):
                divGOmaxphosphodiv.append(par.maxorfphosphodiv)
            divGOminseqdiv.append(par.minorfseqdiv)
            divGOminsytdiv.append(par.minorfsytdiv)
            if hasattr(par, 'minorfphosphodiv'):
                divGOminphosphodiv.append(par.minorfphosphodiv)
        else:
            consGOmaxseqdiv.append(par.maxorfseqdiv)
            consGOmaxsytdiv.append(par.maxorfsytdiv)
            if hasattr(par, 'maxorfphosphodiv'):
                consGOmaxphosphodiv.append(par.maxorfphosphodiv)
            consGOminseqdiv.append(par.minorfseqdiv)
            consGOminsytdiv.append(par.minorfsytdiv)
            if hasattr(par, 'minorfphosphodiv'):
                consGOminphosphodiv.append(par.minorfphosphodiv)





def calGOorf(geneset, orfs, Pars):
    for p in geneset:
        for i,t in enumerate(types):
            if hasattr(Pars[p], 'alignment'):
                if getattr(Pars[p], typesinPars[i]) == str(1):
                    caldivGOorf(geneset, orfs, Pars, t)
                else:
                    calconsGOorf(geneset, orfs, Pars, t)

                    
def caldivGOorf(geneset, orfs, Pars, t, divGOmaxseqdiv = {}, divGOmaxphosphodiv = {}, divGOmaxsytdiv = {}, 
        divGOminseqdiv = {}, divGOminphosphodiv ={}, divGOminsytdiv = {}): 
    if t not in divGOmaxseqdiv.keys():
        divGOmaxseqdiv[t] = []
    divGOmaxseqdiv[t].append(Pars[p].maxorfseqdiv)
    if t not in divGOmaxphosphodiv.keys():
        divGOmaxphosphodiv[t] = []
    if hasattr(Pars[p], 'maxorfphosphodiv'):
        divGOmaxphosphodiv[t].append(Pars[p].maxorfphosphodiv)
    if t not in divGOmaxsytdiv.keys():
        divGOmaxsytdiv[t] = []
    divGOmaxsytdiv[t].append(Pars[p].maxorfsytdiv)
    if t not in divGOminseqdiv.keys():
        divGOmaxseqdiv[t] = []
    divGOminseqdiv[t].append(Pars[p].minorfseqdiv)
    if t not in divGOminphosphodiv.keys():
        divGOminphosphodiv[t] = []
    if hasattr(Pars[p], 'minorfphosphodiv'):    
        divGOminphosphodiv[t].append(Pars[p].minorfphosphodiv)
    if t not in divGOminsytdiv.keys():
        divGOminsytdiv[t] = []
    divGOminsytdiv[t].append(Pars[p].minorfsytdiv)
    printresults(divGOmaxseqdiv, divGOmaxphosphodiv, divGOmaxsytdiv, 
                    divGOminseqdiv, divGOminphosphodiv, divGOminsytdiv)

def calconsGOorf(geneset, orfs, Pars, t, consGOmaxseqdiv = {}, consGOmaxphosphodiv = {},
            consGOmaxsytdiv = {}, consGOminseqdiv = {}, consGOminphosphodiv = {}, consGOminsytdiv = {}):
    if t not in consGOmaxseqdiv.keys():
        consGOmaxseqdiv[t] = []
    consGOmaxseqdiv[t].append(Pars[p].maxorfseqdiv)
    if t not in consGOmaxphosphodiv.keys():
        consGOmaxphosphodiv[t] = []
    if hasattr(Pars[p], 'maxorfphosphodiv'):
        consGOmaxphosphodiv[t].append(Pars[p].maxorfphosphodiv)
    if t not in consGOmaxsytdiv.keys():
        consGOmaxsytdiv[t] = []
    consGOmaxsytdiv[t].append(Pars[p].maxorfsytdiv)
    if t not in consGOminseqdiv.keys():
        consGOminseqdiv[t] = []
    consGOminseqdiv[t].append(Pars[p].minorfseqdiv)
    if t not in consGOminphosphodiv.keys():
        consGOminphosphodiv[t] = []
    if hasattr(Pars[p], 'minorfphosphodiv'):    
        consGOminphosphodiv[t].append(Pars[p].minorfphosphodiv)
    if t not in consGOminsytdiv.keys():
        consGOminsytdiv[t] = []
    consGOminsytdiv[t].append(Pars[p].minorfsytdiv)
    printresults(consGOmaxseqdiv, consGOmaxphosphodiv, consGOmaxsytdiv, 
                    consGOminseqdiv, consGOminphosphodiv, consGOminsytdiv)


def calcdivGO(geneset, orfs, Pars, divGOseqdiv = {}, divGOphosphodiv = {}, divGOsytdiv = {}, consGOseqdiv = {}, 
        consGOphosphodiv = {}, consGOsytdiv = {}):
    for p in geneset:
        for i, t in enumerate(types):
            if hasattr(Pars[p], 'alignment'):
                if getattr(Pars[p], typesinPars[i]) == str(1):
                    if t not in divGOseqdiv.keys():
                        divGOseqdiv[t] = []
                    divGOseqdiv[t].append(Pars[p].seqdiff)
                    if t not in divGOphosphodiv.keys():
                        divGOphosphodiv[t] = []
                    divGOphosphodiv[t].append(1-Pars[p].phosphositecons)
                    if t not in divGOsytdiv.keys():
                        divGOsytdiv[t] = []
                    divGOsytdiv[t].append(Pars[p].sytdiff)
                else:
                    if t not in consGOseqdiv.keys():
                        consGOseqdiv[t] = []
                    consGOseqdiv[t].append(Pars[p].seqdiff)
                    if t not in consGOphosphodiv.keys():
                        consGOphosphodiv[t] = []
                    consGOphosphodiv[t].append(1-Pars[p].phosphositecons)
                    if t not in consGOsytdiv.keys():
                        consGOsytdiv[t] = []
                    consGOsytdiv[t].append(Pars[p].sytdiff)
    printresults(divGOseqdiv, divGOphosphodiv, 
            divGOsytdiv, consGOseqdiv, consGOphosphodiv, consGOsytdiv)

def printresults(divGOseqdiv, divGOphosphodiv, divGOsytdiv, consGOseqdiv,
                      consGOphosphodiv, consGOsytdiv):
    f = fisher.FisherExactTest()
    for t in types:
        means = [numpy.mean(divGOseqdiv[t]), numpy.mean(divGOsytdiv[t]), numpy.mean(consGOseqdiv[t]), numpy.mean(consGOsytdiv[t])]
        meansphosphos= [numpy.mean(divGOseqdiv[t]), numpy.mean(divGOphosphodiv[t]), numpy.mean(consGOseqdiv[t]), numpy.mean(consGOphosphodiv[t])]
        print "means:\n"
        print "divGOseqdiv, divGOsytdiv, consGOseqdiv, consGOsytdiv"
        print means
        f.print_report(means[0], means[1], means[2], means[3])
        "means_phosphos:\n"
        print "divGOseqdiv, divGOphosphodiv, consGOseqdiv, consgophosphodiv"
        print meansphosphos
        f.print_report(meansphosphos[0], meansphosphos[1], meansphosphos[2], meansphosphos[3])







                    
                

