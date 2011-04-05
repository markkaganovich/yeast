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


def calcdivGO(geneset, divGOseqdiv, divGOphosphodiv, divGOsytdiv, consGOseqdiv, 
        consGOphosphodiv, consGOsytdiv, orfs, Pars):
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







                    
                

