'''
yeast gene duplication results in new function, bitches

the existence of orfs and Pars is assumed
'''
import fisher
import numpy
import operator

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

maxseqdiv = []
minseqdiv = []
maxphosphodiv = []
minphosphodiv = []


for p in geneset:
    par = Pars[p]
    if hasattr(par, 'maxorfseqdiv'):
        maxseqdiv.append(par.maxorfseqdiv)
        if hasattr(par, 'maxorfphosphodiv'):
            maxphosphodiv.append(par.maxorfphosphodiv)
        minseqdiv.append(par.minorfseqdiv)
        if hasattr(par, 'minorfphosphodiv'):
            minphosphodiv.append(par.minorfphosphodiv)

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


for p in geneset:
    par = Pars[p]
    if hasattr(par, 'maxorfseqdiv'):
        maxseqdiv.append(par.maxorfseqdiv)
        if hasattr(par, 'maxorfphosphodiv'):
            divGOmaxphosphodiv.append(par.maxorfphosphodiv)
        minseqdiv.append(par.minorfseqdiv)
        if hasattr(par, 'minorfphosphodiv'):
            minphosphodiv.append(par.minorfphosphodiv)


oneisaTF = []
bothareTFs = []
for p in Pars.keys():
    par = Pars[p]
    if par.orfs[0] in orfs.keys() and par.orfs[1] in orfs.keys():
        o1 = orfs[par.orfs[0]]
        o2 = orfs[par.orfs[1]]
        if sum([o1.TF, o2.TF]) == 2:
            bothareTFs.append(p)
        if sum([o1.TF, o2.TF]) == 1:
            oneisaTF.append(p)
    
TFparsfromWGD =[]
TFparsfromSSD = []
for p in bothareTFs:
    if Pars[p].event == 'WGD':
        TFparsfromWGD.append(p)
    else:
        TFparsfromSSD.append(p)

oneisaTFfromWGD = []
oneisaTFfromSSD = []
for p in oneisaTF:
    if Pars[p].event == 'WGD':
        oneisaTFfromWGD.append(p)
    else:
        oneisaTFfromSSD.append(p)

phosphosites = []
for o in TFs:
    if o not in orfs.keys():
        o
        continue
    if hasattr(orfs[o], 'phosphosite'):
        phosphosites.append(len(orfs[o].phosphosite))

phosphosites = sum(phosphosites)

simlength = 1000
simphosphosums = []

for j in range(0,simlength):
    phosphosites = []
    for i in range(0,len(TFs)):
        n = numpy.random.randint(len(orfs.keys()))
        orf = orfs[orfs.keys()[n]]
        if hasattr(orf, 'phosphosite'):
            phosphosites.append(len(orf.phosphosite))

    simphosphosums.append(sum(phosphosites))

# make a geneset of just TFs that have phosphosites
geneset = []
allTFpars = oneisaTF + bothareTFs
for p in allTFpars:
    o1 = orfs[Pars[p].orfs[0]]
    o2 = orfs[Pars[p].orfs[1]]
    if hasattr(o1, 'phosphosite') and hasattr(o2, 'phosphosite'):
        geneset.append(p)


# Rick Young 2002 Science TF binding regions
youngTFslist #second line of file
TFnames #Jen's TF list converted to gene names
file = open('./data/binding_by_gene.tsv')
lines = file.readlines()
file.close()
TFbinding ={}
for tf in TFs:
    TFbinding[tf] = []

for line in lines:
    l = line.strip('\n').strip('\r').split('\t')
    gene = l[0]
    for i in range(4, len(l)):
        if l[i] !='' and float(l[i]) < .05:
            youngTF = youngTFslist[i]
            if youngTF in TFnames:
                TFbinding[TFs[TFnames.index(youngTF)]].append(gene)

indofp = [i for (i,j) in sorted(enumerate(phosphodiff), key=operator.itemgetter(1), reverse = True)]














                
                
                
                
                
                
                
                
                
                
for p in geneset:
    par = Pars[p]
    if hasattr(par, 'maxorfseqdiv') and hasattr(par, 'minorfphosphodiv'):
        here.append(p)
        if getattr(par, type) == str(1):
            divGOmaxseqdiv.append(par.maxorfseqdiv)
            divGOmaxsytdiv.append(par.maxorfsytdiv)
            divGOmaxphosphodiv.append(par.maxorfphosphodiv)
            divGOminseqdiv.append(par.minorfseqdiv)
            divGOminsytdiv.append(par.minorfsytdiv)
            divGOminphosphodiv.append(par.minorfphosphodiv)
        else:
            consGOmaxseqdiv.append(par.maxorfseqdiv)
            consGOmaxsytdiv.append(par.maxorfsytdiv)
            consGOmaxphosphodiv.append(par.maxorfphosphodiv)
            consGOminseqdiv.append(par.minorfseqdiv)
            consGOminsytdiv.append(par.minorfsytdiv)
            consGOminphosphodiv.append(par.minorfphosphodiv)




    for p in geneset:
        for i,t in enumerate(types):
            if hasattr(Pars[p], 'alignment'):
                if getattr(Pars[p], typesinPars[i]) == str(1):
                    caldivGOorf(geneset, orfs, Pars, t)
                else:
                    calconsGOorf(geneset, orfs, Pars, t)

                    
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







                    
                

