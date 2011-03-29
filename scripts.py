from numpy import *
from helpers import *
from YeastPhosphoProject import *
from pylab import *

species = LoadSpecies()
pre_wgd = [6, 7, 8, 10]
post_wgd = [1,2,3,4,5,9]

class Protein:
    def __init__(self, index, alignment, phosphosites, scergenes):
                self.index = index

                self.name = scergenes[index]
                self.prtn_len = len(filter(lambda x: x!= '-', alignment[index][0]))
                self.num_phosphos = len(phosphosites[index])
            
                self.cons_sites_in_ortho = sum(map(lambda x: conservation(index, x, [8], alignment), phosphosites[index]))
                self.cons_aas_in_ortho = all_cons(index, alignment, [0,8])
                
                
                [self.numsty, self.conssty] = parSTYcons(alignment[index], 8)



class ParalogPairs:
    def __init__(self, index, par_alignment, par_phosphosites, alignment, phosphosites, wgd=1):
        import simplejson as json
        from scripts import all 

        genenamesfile = open('genenamesfile')
        genenames = json.load(genenamesfile)
        genenamesfile.close()

        ##go categories
        file = open('goComponent')
        goComponent = json.load(file)
        file.close()
        file = open('goFunction')
        goFunction = json.load(file)
        file.close()
        file = open('goProcess')
        goProcess = json.load(file)
        file.close()
        file = open('goComponent_rev')
        goComponent_rev= json.load(file)
        file.close()
        file = open('goProcess_rev')
        goProcess_rev = json.load(file)
        file.close()
        file = open('goFunction_rev')
        goFunction_rev = json.load(file)
        file.close()
        #################################


        self.wgd = wgd    
        print index
        self.index = index
        index1 = index[0]
        index2 = index[1]
        self.name = [genenames[index1], genenames[index2]]
        self.prtn_len = map(lambda x: double(len(filter(lambda x: x!= '#', par_alignment[genenames[index1]+'#'+genenames[index2]][x]))), [0, 1])
        
        name1 = genenames[index1]
        name2 = genenames[index2]
        name = name1+'#'+name2
        self.phosphoevents = [par_phosphosites[name1+'#'+name2][0], par_phosphosites[name1+'#'+name2][1]]
        self.num_phosphoevents = len(self.phosphoevents[0])+len(self.phosphoevents[1])
        self.events_in_common = list(set(par_phosphosites[name1+'#'+name2][0]) & set(par_phosphosites[name1+'#'+name2][1]))
        self.events_differ = filter(lambda x: x not in self.events_in_common, self.phosphoevents[0]+self.phosphoevents[1])
    
        self.phosphositesInPara = [[],[]]
        for i in par_phosphosites[name][0]:
            self.phosphositesInPara[1].append(par_alignment[name][1][i])
        for i in par_phosphosites[name][1]:
            self.phosphositesInPara[0].append(par_alignment[name][0][i])

        # paralogs in same GO slim categories?:
        if goComponent.get(name1) != 'cellular_component' and goComponent.get(name2) != 'cellular_component':
            if goComponent.get(name1) == goComponent.get(name2):
                self.sameComponent = True
            else:
                self.sameComponent = False
        else:
             self.sameComponent = 'Undefined'

            
        if goProcess.get(name1) != 'biological_process' and goProcess.get(name2) != 'biological_process':
            if goProcess.get(name1) == goProcess.get(name2):
                self.sameProcess = True
            else:
                self.sameProcess = False
        else:
             self.sameProcess = 'Undefined'


        if goFunction.get(name1) != 'molecular_function' and goFunction.get(name2) != 'molecular_function':
            if goFunction.get(name1) == goFunction.get(name2):
                self.sameFunction = True
            else:
                self.sameFunction = False
        else:
             self.sameFunction = 'Undefined'

        self.cons_sites_in_ortho1 = sum(map(lambda x: conservation(index1, x, [8], alignment), phosphosites[index1]))
        self.cons_sites_in_ortho2 = sum(map(lambda x: conservation(index2, x, [8], alignment), phosphosites[index2]))
        self.cons_seq1 = sum(map(lambda x: conservation(index1, x, [8], alignment), range(0, len(alignment[index1][0]))))/double(self.prtn_len[0])
        self.cons_seq2 = sum(map(lambda x: conservation(index2, x, [8], alignment), range(0, len(alignment[index2][0]))))/double(self.prtn_len[1])
        self.ortho_sites_in_common = len(list(set(phosphosites[index1]) & set(phosphosites[index2])))
        self.num_ortho_phosphoevents = len(phosphosites[index1]) + len(phosphosites[index2])
    
        self.cons_sites_ortho_norm = (self.cons_sites_in_ortho1 + self.cons_sites_in_ortho2) / self.num_ortho_phosphoevents

        self.num_conserved_sites_bysequence1 = sum(map(lambda x: conservation(index1, x, [1], par_alignment[name1+'#'+name2]), par_phosphosites[name1+'#'+name2][0]))
        self.num_conserved_sites_bysequence2 = sum(map(lambda x: conservation(index1, x, [1], par_alignment[name1+'#'+name2]), par_phosphosites[name1+'#'+name2][1]))
        self.total_cons_phosphos = self.num_conserved_sites_bysequence1 + self.num_conserved_sites_bysequence2 - self.events_in_common
        if sum(self.num_phosphoevents) == 0:
            self.total_cons_phosphos_norm = -.5
        else:
            self.total_cons_phosphos_norm = self.total_cons_phosphos/sum(self.num_phosphoevents)
        self.seq_cons = [all_cons_pars(index1, par_alignment[name1+'#'+name2], [1]), all_cons_pars(index2, par_alignment[name1+'#'+name2], [1])]


        [self.numsty, self.conssty] = parSTYcons(par_alignment[name])
        
        [self.prtn_len, self.cons_aa] = all_cons_par2(par_alignment[name])       
   
def parsrun(parsls, DICparsalignment, parphosphosites_dic, alignment, phosphosites):
    return map(lambda y:  map(lambda x: ParalogPairs(x, DICparsalignment, parphosphosites_dic, alignment, phosphosites, 1), y), parsls)
        


class speciesParalogPairs:
    def __init__(self, sp, orthos, index1, index2, par_alignment, par_phosphosites, alignment, phosphosites, wgd=1):
        import simplejson as json
        from scripts import all 

        genenames = LoadfromFile('../genomes/'+sp+'Genes')

        self.wgd = wgd    

        species = ['Spar', 'Smik', 'Skud', 'Sbay', 'Scas', 'Sklu', 'Klac', 'Calb', 'Cgla',  'Spom']
        
        speciesindex = species.index(sp)
        orthoindex1 = orthos[index1][speciesindex]
        orthoindex2 = orthos[index2][speciesindex]

        self.index = [index1, index2]
        self.name = [genenames[orthoindex1], genenames[orthoindex2]]
#self.prtn_len = map(lambda x: double(len(filter(lambda x: x!= '-', par_alignment[genenames[orthoindex1]+'#'+genenames[orthoindex2]][x]))), [0, 1])
        
        name1 = genenames[orthoindex1]
        name2 = genenames[orthoindex2]
        if name1+'#'+name2 in par_phosphosites.keys():
            self.phosphoevents =par_phosphosites[name1+'#'+name2]
            print self.phosphoevents
            self.num_phosphoevents = len(self.phosphoevents)
             
            self.num_conserved_sites_bysequence = sum(map(lambda x: conservation(index1, x, [1], par_alignment[name1+'#'+name2]), par_phosphosites[name1+'#'+name2]))
            if sum(self.num_phosphoevents) == 0:
                self.total_cons_phosphos_norm = -.5
            else:
                self.total_cons_phosphos_norm = self.total_cons_sites_bysequence/sum(self.num_phosphoevents)
        self.cons_sites_in_ortho1 = sum(map(lambda x: conservation(index1, x, [speciesindex+1], alignment), phosphosites[index1]))
        self.cons_sites_in_ortho2 = sum(map(lambda x: conservation(index2, x, [speciesindex+1], alignment), phosphosites[index2]))
        self.cons_seq1 = sum(map(lambda x: conservation(index1, x, [speciesindex+1], alignment), range(0, len(alignment[index1][0]))))/double(self.prtn_len[0])
        self.cons_seq2 = sum(map(lambda x: conservation(index2, x, [speciesindex+1], alignment), range(0, len(alignment[index2][0]))))/double(self.prtn_len[1])
        self.ortho_sites_in_common = len(list(set(phosphosites[index1]) & set(phosphosites[index2])))
        self.num_ortho_phosphoevents = len(phosphosites[index1]) + len(phosphosites[index2])
    
        self.cons_sites_ortho_norm = (self.cons_sites_in_ortho1 + self.cons_sites_in_ortho2) / self.num_ortho_phosphoevents

#self.seq_cons = [all_cons_pars(index1, par_alignment[name1+'#'+name2], [1]), all_cons_pars(index2, par_alignment[name1+'#'+name2], [1])]
        [self.prtn_len, self.cons_aa] = all_cons_pars2(par_alignment[name])       

def temp(l_high, cons_phosphos_sorted_f, scergenes):
    l_high1 = []
    l_high2= []
    for i in range(0,len(l_high)):
        if cons_phosphos_sorted_f[106:151][i].num_phosphoevents[0] > cons_phosphos_sorted_f[106:151][i].num_phosphoevents[1]:
            l_high1.append(scergenes[l_high[i][0]])
            l_high2.append(scergenes[l_high[i][1]])
        else: 
            l_high1.append(scergenes[l_high[i][1]])
            l_high2.append(scergenes[l_high[i][0]])

    return [l_high1, l_high2]

def count_pars(dir='../genomes/WGD_pars_alignments'):
    files = os.listdir(dir)
    dic={}
    for file in files:
        genes=[]
    
        genes.append(file.split('_')[0])
        genes.append(file.split('_')[1].split('.')[0])
        for gene in genes:
            if gene not in dic:
                dic.update({gene : 1})
            else:
                dic.update({gene : dic[gene]+1})

    return dic
    
    
def calc_sites_by_protein(phosphosites, alignment):
    protein_sites = []
    for i in range(0,len(phosphosites)):
        if alignment[i][0] != []:
            prtn_len = double(len(filter(lambda x: x != '-', alignment[i][0])))
            protein_sites.append(Protein(i, prtn_len, len(phosphosites[i]), len(filter(lambda x: x == 'S' or x=='T', alignment[i][0]))/prtn_len))

    return protein_sites  

def pick_random_genes(n):
    randgenes=[]
    for i in range(0,n):
        randgenes.append(int(rand()*5800))

    return randgenes

def conservation(gene, index, species_set, alignment):
    if len(alignment) > 2:
        if alignment[gene][0] == [] or alignment[gene][1] == []:
            return -1
        l =filter(lambda x: alignment[gene][x] !=[] and alignment[gene][x][index] == alignment[gene][0][index]  , species_set)
        cr = double(len(l))/len(species_set)
        return cr
    else:
        if alignment[0] == [] or alignment[1] == []:
            return -1
        l = filter(lambda x: alignment[x] !=[] and alignment[x][index] == alignment[0][index], species_set)
        cr = double(len(l))/len(species_set)
        return cr

def avg_conservation_per_prtn(genes, phosphosites, species_set, alignment):
    avg_phosphosites=[]
    complete_alignment_genes = list(set(complete_aligns(alignment, species_set)) & set(genes))
    for i in complete_alignment_genes:
        print i
        avg_phosphosites.append(average(map(lambda x: conservation(i, x, species_set, alignment), phosphosites[i])))
    avg = filter(lambda x: str(x)!= 'nan', avg_phosphosites)
    
    return avg
        

def complete_aligns(alignment, species_set):
    complete_aligns=[]
    for i in range(0, len(alignment)):
        complete = filter(lambda x: alignment[i][x] != [], species_set)
        if len(complete) == len(species_set):
            complete_aligns.append(i)
    return complete_aligns

def paralog_pair_cons(pair1, pair2, species_set, alignment, phosphosites):
    new_pair1=[]
    new_pair2=[]
    diffs=[]
    complete_align = complete_aligns(alignment, species_set)
    for i in range(0, len(pair1)):
        t1= all_cons(pair1[i], alignment, species_set)
        t2 = all_cons(pair2[i], alignment, species_set)
        a1 = average(map(lambda x: conservation(pair1[i], x, species_set, alignment), phosphosites[pair1[i]]))
        a2 = average(map(lambda x: conservation(pair2[i],x,species_set, alignment), phosphosites[pair2[i]]))
        if a1/t1 <= a2/t2:    
            new_pair1.append(pair1[i])
            new_pair2.append(pair2[i])
            diffs.append(a2/t2-a1/t1)
        else:
            new_pair1.append(pair2[i])
            new_pair2.append(pair1[i])
            diffs.append(a1/t1-a2/t2)
    return [new_pair1, new_pair2, diffs]

def all_cons(gene, alignment, species_set):
    
    if alignment[gene] == []:
        return 0
    return len([j for j in range(0, len(alignment[gene][0])) if len(filter(lambda x: alignment[gene][x] != [] and alignment[gene][x][j] != '-' and alignment[gene][x][j] == alignment[gene][0][j],species_set)) == len(species_set)-1])


def all_cons_pars(gene, alignment, par):
    if len(alignment) > 2: 
        if alignment[gene] == []:
            return 0
        return len([j for j in range(0, len(alignment[gene][0])) if len(filter(lambda x: alignment[gene][x] != [] and alignment[gene][x][j] != '-' and alignment[gene][x][j] == alignment[gene][0][j],par)) == 1])/ double(len(filter(lambda x: x!='-', alignment[gene][0])))
    else:
        return len([j for j in range(0, len(alignment[0])) if len(filter(lambda x: alignment[x] != [] and alignment[x][j] != '-' and alignment[x][j] == alignment[0][j],par)) == 1])/ double(len(filter(lambda x: x!='-', alignment[0])))

# a different accounting        
def all_cons_par2( alignment):
    if alignment != [] and alignment[0] == [] or alignment[1] == []:
        prtn_len = 0
        cons_aa = 0
    else:
        prtn_len = len(alignment[0])
        cons_aa = len([j for j in range(0, len(alignment[0])) if alignment[0][j] == alignment[1][j]])

    return [prtn_len, cons_aa]


        

def parSTYcons(alignedseq, species = 1):
    numSTY=0
    consSTY=0

    if len(alignedseq) >= 2:
        if alignedseq[0] ==[] or alignedseq[species] == []:
            return [0,0]
        stypos = list(set([j for j in range(0,len(alignedseq[0])) if alignedseq[0][j] == 'S' or alignedseq[0][j] == 'T' or alignedseq[0][j] == 'Y'] + [j for j in range(0,len(alignedseq[species])) if alignedseq[species][j] == 'S' or alignedseq[species][j] == 'T' or alignedseq[species][j] == 'Y']))
        consSTY = len(filter(lambda x: alignedseq[0][x] == alignedseq[species][x], stypos))

        numSTY = len(stypos)

    return [numSTY, consSTY]
def parphosphopos(wgd, phosphosites, alignment, species, orthos):
  
    parphosphosites_dic = {}
    specieslist = ['Spar', 'Smik', 'Skud', 'Sbay', 'Scas', 'Sklu', 'Klac', 'Calb', 'Cgla',  'Spom']
    spindex = int(specieslist.index(species))
    print spindex
    genenames = LoadfromFile('../genomes/'+species+'Genes')
    file = open('./'+species+'_parsalignment')
    parsalign = json.load(file) 
    file.close()

    for x in wgd : 
        if phosphosites[x[0]] != [] and orthos[x[0]][spindex] != [] and orthos[x[1]][spindex] != []:
            name = str(genenames[orthos[x[0]][spindex]])+'#'+str(genenames[orthos[x[1]][spindex]])
            if name not in parsalign.keys():
                print False
                print name
            else:
                temp =[]
                for p in phosphosites[x[0]]:
                    temp.append(TranslatePos2Align(Align2Pos(p,alignment, 0,1), parsalign[name][0]))
                parphosphosites_dic.update({name : temp})

    return parphosphosites_dic
def log_regression_paralogs(geneset, wgd, alignment, phosphosites, species = 8):
    conservation = []
    pho =[]
    inwgd =[]
    numgenes=0
    for gene in geneset:
        print gene
        if alignment[gene][0] != [] and alignment[gene][species]!= []:
            numgenes = numgenes+1
            for j in range(0, len(alignment[gene][0])):
                if not (alignment[gene][0][j] == '-' and alignment[gene][species][j]=='-'):
                    if alignment[gene][0][j] == alignment[gene][species][j]:
                        conservation.append(1)
                    else:
                        conservation.append(0)
                    if j in phosphosites[gene]:
                        pho.append(1)
                    else:
                        pho.append(0)
                    if gene in wgd:
                        inwgd.append(1)
                    else:
                        inwgd.append(0)
    print numgenes
    return [conservation, pho, inwgd]
#----------------------------------------------------------------------------
# stuff for paralog comparisons pair:pair
# ---------------------------------------------------------------------------
def log_regression_paralogs_pair(f,  alignment, phosphosites):
    conservation = []
    pho =[]
    inwgd =[]
    numgenes=0
    for pair in f:
        name = str(pair.name[0])+'_'+str(pair.name[1])
        print name
        if alignment[name][0] != [] and alignment[name][1]!= []:
            numgenes = numgenes+1
            for j in range(0, len(alignment[name][0])):
                if not (alignment[name][0][j] == '-' and alignment[name][1][j]=='-'):
                    if alignment[name][0][j] == alignment[name][1][j]:
                        conservation.append(1)
                    else:
                        conservation.append(0)
                    if j in phosphosites[name][0]+phosphosites[name][1]:
                        pho.append(1)
                    else:
                        pho.append(0)
                    if pair.wgd == 1:
                        inwgd.append(1)
                    else:
                        inwgd.append(0)
    print numgenes
    return [conservation, pho, inwgd]


def log_regression_paralogs_pair_truephosphosite(f_pars, alignment, phosphosites, windowsize = 5):
    eventconservation = [] # whether the actual phosphosite event is conserved... 
    conservationinwindow = []  # how many amino acids are conserved in the window
    inwgd = [] # binary, in wgd or not
    siteconservation =[] # binary, is the seq  site conserved or not

    for pair in f_pars:
        name = str(pair.name[0])+'#'+str(pair.name[1])
        print name
        window = windowsize
        if alignment[name][0] != [] and alignment[name][1] != []:
            for j in list(set(pair.phosphoevents[0] +pair.phosphoevents[1])):
                if j in pair.events_in_common:
                    eventconservation.append(1)
                else:
                    eventconservation.append(0)
                window = windowsize
                tempcons = 0
                for k in range(1, window):   # going right
                    if j+k>=0 and j+k<min(len(alignment[name][0]), len(alignment[name][1])):
                        if alignment[name][0][j+k] == '-' and alignment[name][1][j+k] == '-':
                            window = window +1
                        else:
                            if alignment[name][0][j+k] == alignment[name][1][j+k]:
                                tempcons = tempcons+1
                window = windowsize
                for k in range(1, window): # going to left 
                    if j-k >= 0 :
                        if alignment[name][0][j-k] == '-' and alignment[name][1][j-k] == '-':
                            window = window+1
                        else:
                            if alignment[name][0][j-k] == alignment[name][1][j-k]:
                                tempcons = tempcons+1
                conservationinwindow.append(tempcons)

                if alignment[name][0][j] == alignment[name][1][j]:
                    siteconservation.append(1)
                else:
                    siteconservation.append(0)
                
                if pair.wgd == 1:
                    inwgd.append(1)
                else:
                    inwgd.append(0)


            

    return [eventconservation, conservationinwindow, siteconservation, inwgd]
                        

    


def go_analysis(new_pair1, new_pair2, genes, GOref, GOslim):
    
    same = 0
    different = 0
    for i in genes:
        if GOslim[new_pair1[i]][2] == GOslim[new_pair2[i]][2]: 
            same = same+1
        else:
            different = different+1

    return [different, same]
