# April 29, 2009
#
# Mark Kaganovich
#
# Yeast phosphorylation site conservation project scripts
#
import sys
import os
import shutil
import simplejson as json
from helpers import *

# Global variables
SCer = 0

NUM_COMPARISONS = 44
CONS = 0
SYT = 1
ED = 2
DIV = 3
PW_1 = 4            # conserved completely in pw comparison
PW_2 = 5
PW_3 = 6
PW_4 = 7
PW_5 = 8
PW_6 = 9
PW_7 = 10
PW_8 = 11
PW_9 = 12
PW_10 = 13
PW_D_1 = 14         # not completely conserved
PW_D_2 = 15
PW_D_3 = 16
PW_D_4 = 17
PW_D_5 = 18
PW_D_6 = 19
PW_D_7 = 20
PW_D_8 = 21
PW_D_9 = 22
PW_D_10 = 23
PW_SYT_1 = 24       # not conserved, but mutated within SYT
PW_SYT_2 = 25
PW_SYT_3 = 26
PW_SYT_4 = 27
PW_SYT_5 = 28
PW_SYT_6 = 29
PW_SYT_7 = 30
PW_SYT_8 = 31
PW_SYT_9 = 32
PW_SYT_10 = 33
PW_SYT_D_1 = 34     # not conserved AND not mutated with SYT
PW_SYT_D_2 = 35
PW_SYT_D_3 = 36
PW_SYT_D_4 = 37
PW_SYT_D_5 = 38
PW_SYT_D_6 = 39
PW_SYT_D_7 = 40
PW_SYT_D_8 = 41
PW_SYT_D_9 = 42
PW_SYT_D_10 = 43


F = 0
P = 1
C = 2


# Alignment(GENE, SPECIES)

#parse out only WGD genes from 'proteins' folder downloaded from Kellis et al paper on duplication
def OnlyWGD(scergenes):
    files = os.listdir('../genomes/protein')
    WGDpairs=[]
    WGDind = []
    wgd_after1 = []
    wgd_after2 =[]
    print(len(files))
    for i in range(0, len(files)):
        names = files[i].split('_')
        WGDpairs.append(names[0])
        WGDpairs.append(names[1][0:-5])
        wgd1_temp = FindinList(names[0], scergenes)
        wgd2_temp = FindinList(names[1][0:-5], scergenes)
        if wgd1_temp != -1:
            wgd_after1.append(FindinList(names[0], scergenes))
            wgd_after2.append(FindinList(names[1][0:-5], scergenes))
    return [wgd_after1, wgd_after2]


# intersect output of onlywgd with paralog_genes        
def after_wgd_pairs(pair1, pair2, paralog_genes):
    result_pair1=[]
    result_pair2=[]
    for i in range(0, len(pair1)):
        if pair1[i] in paralog_genes and pair2[i] in paralog_genes:
            result_pair1.append(pair1[i])
            result_pair2.append(pair2[i])
#return filter(lambda x,y: x in paralog_genes and y in paralog_genes, pair1 pair2)
    return [result_pair1, result_pair2]


        
def onlywgd2(pars, wgdind):
    onlywgd = []
    otherpars = []
    for i in range(0, len(pars)):
        if FindinListIntegers(pars[i], wgdind) != -1:
            onlywgd.append(pars[i])
        else:
            otherpars.append(pars[i])

    return [onlywgd, otherpars]

    
# new align parse: July 14, 2009. Parse outputs that have entire seq in one line
def AlignParse(dir):
    files = os.listdir('../genomes/'+dir)
    Alignment = []

    species = GenSpecies()
    ScerGenes = LoadfromFile('../genomes/ScerGenes')
        
    for i in range(0, len(ScerGenes)):
        Alignment.append([])
        for j in range(0, len(species)):
            Alignment[i].append([])

    for i in range(1,len(files)):
        name = files[i].split('.')
        name = name[0]
        print name
        index = FindinList(name, ScerGenes)
            
        File=open('../genomes/'+dir+'/'+files[i], 'r')
        line = 1
        seq = []
        line = File.readline()
        while (line):
            if line[0] == '>':
                if seq != []:
                    if t != -1 and line != '' and index != -1:
                        Alignment[index][t].extend(seq)
                    seq = []
                sp = str(line.split('\t')[0])[1:5]
                t = FindinList(sp, species)
            else: 
                seq.extend(line.split('\n')[0])
            line = File.readline()
        
        Alignment[index][t].extend(seq)

    return Alignment


# parse paralog-paralog alignments, in parsalignments folder
def ParParse(ScerGenes, dir ):
    files = os.listdir('../genomes/'+dir)
    Alignment = {}

    print files
    ScerGenes = json.load(open('genenamesfile'))    
    
#for i in range(0, len(ScerGenes)):
#       Alignment.append([])
#       for j in range(0, 2):
#           Alignment[i].append([])

    for i in range(1,len(files)):
        name = files[i].split('.')
        name = name[0]
#name = files[i][0:-7]
        name1 = name.split('_')[0]
        name2 = name.split('_')[1]
#name1 = 'Scas'+name.split('Scas')[1][0:-1]
#       name2 = 'Scas'+name.split('Scas')[2]
        print name1
        print name2
        index1 = FindinList(name1, ScerGenes)
        index2 = FindinList(name2, ScerGenes)
            
        File=open('../genomes/'+dir+'/'+files[i], 'r')
        line = 1
        seq = [[],[]]
        line = File.readline()
        t = -1
        while (line):
            if line[0] == '>':
                t=t+1    
            else: 
                seq[t].extend(line.split('\n')[0])
            
            line = File.readline()
        
        if seq[0] == [] or seq[1] == []:
            Alignment.update({name1+'#'+name2:[[],[]]})
        else:
            Alignment.update({name1+'#'+name2:seq})
        
    return Alignment


def checkalignment(alignment):

    species = GenSpecies()
    full = [0]*len(species)

    for i in range(0, len(alignment)):
        for s in range(0, len(species)):
            if alignment[i][s] != [] and len(alignment[i][s])>10:
                full[s] = full[s]+1

    return full


#Input: tab delimited files with phosphorylation site position in S. Cer for every gene
# filename is: './PhosphorylationSites.tab
def ParsePhosphoSiteFiles(filename, genes, alignment):
    
    SCerPhosphoSites = []
    for i in range(0, len(genes)):
        SCerPhosphoSites.append([])

    File = open('./' + filename, 'r')
    line = File.readline()
    while (line):
        line = File.readline()
        geneindex = FindinList(str(line.split('\t')[0]).upper(),genes)
        if geneindex == -1:
            print str(line.split('\t')[0])
        
        if geneindex != -1 and line != '':
            SCerPhosphoSites[geneindex].append(TranslatePos2Align(int(line.split('\t')[1]), alignment[geneindex][0]))
            print int(line.split('\t')[1])
#           print TranslatePos2Align(int(line.split('\t')[1]), alignment[geneindex][0])
        else:
            print line.split('\t')[0]

    return SCerPhosphoSites


def ParseParPhosphoSiteFiles(filename, genes, alignment, orthos):
    
    SCerPhosphoSites = []
    for i in range(0, len(genes)):
        SCerPhosphoSites.append([])

    File = open('./' + filename, 'r')

    lines = File.readlines()
    for line in lines:
        gene = str(line.split('\t')[0])
        if gene in genes:
            print gene
            SCerPhosphoSites[genes.index(gene)].append(int(line.split('\t')[1]))

    phosphosites_dic={}
    for key in alignment.keys():
        print key
        phosphosites_dic.update({key:[[],[]]})
        name1 = key.split('#')[0]
        name2 = key.split('#')[1]
        names = [name1, name2]
        for i in range(0, len(names)):
            print names[i]
            for pos in SCerPhosphoSites[genes.index(names[i])]:
                print TranslatePos2Align(pos, alignment[key][i])
                temp=[]
                temp.extend(phosphosites_dic[key])
                temp[i].append(TranslatePos2Align(pos, alignment[key][i]))
                phosphosites_dic.update({key:temp})

    return phosphosites_dic

def FindinList(item, list, num = 1):
    count = 0
    if list != []:
        for i in range(0,len(list)):
            if item in list[i]:
                count = count+1
                if count == num:
                    return i
            else:
                result = -1
    else:
         result = -1
            
    return result

def FindinListIntegers(item, list, num = 1):
    count = 0
    if list != []:
        for i in range(0,len(list)):
            if item ==  list[i]:
                count = count+1
                if count == num:
                    return i
            else:
                result = -1
    else:
         result = -1

    return result

def FindinListofLists(item, listoflists):
    count = 0
    if list != []:
        for i in range(0, len(listoflists)):
            if FindinListIntegers(item, listoflists[i]) != -1:
                return 1

def LoadSpecies():
    species = ['Scer' , 'Spar' , 'Smik', 'Skud', 'Sbay', 'Scas', 'Sklu', 'Klac', 'Calb', 'Cgla', 'Spom']

    return species

#parse from tab delimited file
def BuildGOterms(filename, genes):

    
    File = open('./' + filename, 'r')
    
    GO_slim = []
    for i in range(0, len(genes)):
        GO_slim.append([])
        for j in range(0,3):
            GO_slim[i].append([])

    line = File.readline()
    while (line):
        geneindex = FindinList(str(line.split('\t')[0]).upper(), genes)
        if geneindex != -1 and line != '':
            name = str(line.split('\t')[3])
            cat = str(line.split('\t')[2])
            if cat == 'F':
                GO_slim[geneindex][F].append(name[0:len(name)-2])
            if cat == 'P':
                GO_slim[geneindex][P].append(name[0:len(name)-2])
            if cat == 'C':
                GO_slim[geneindex][C].append(name[0:len(name)-2])
        
        line = File.readline()

    return GO_slim

# build reference where words are indeces
def BuildGOref(GO_slim_names):

    GOref = []    
    for j in range(0,3):
        GOref.append([])

    for i in range(0,len(GO_slim_names)):
        for j in range(0,3):
            for k in range(0,len(GO_slim_names[i][j])):
                name = GO_slim_names[i][j][k]
#print FindinList(name, GOref[j])
                if FindinList(name, GOref[j]) == -1:
                    GOref[j].append(name)

    return GOref
            
def LoadGOref(genes):

    GO = BuildGOterms('go_slim_mapping.txt', genes)

    GOref = BuildGOref(GO)

    return GOref

def LoadGOterms(genes):

    GO = BuildGOterms('go_slim_mapping.txt', genes)
    GOref = BuildGOref(GO)

    GOslim = []
    for i in range(0, len(GO)):
        GOslim.append([])
        for j in range(0,3):
            GOslim[i].append([])
            for k in range(0,len(GO[i][j])):
                GOslim[i][j].append([])
                GOslim[i][j][k] = FindinList(GO[i][j][k], GOref[j])

    return GOslim


#Input: sequence position in a species and the relevant sequence
#Ouput: This position in the alignment (i.e. add up gaps)
def TranslatePos2Align(SeqPos, seq ):
    count=0
    for i in range(0,len(seq)):
        if seq[i] != '-':
            count = count+1
            if count == SeqPos:
                return i
    return -1

def Align2Pos(seqpos, alignment, gene, s):
    count = 0
    for i in range(0, len(alignment[gene][s])):
        if alignment[gene][s][i] != '-':
            count = count+1
            if i == seqpos:
                return count

    return -1
#Input: tab delimited SCer phosphosite file with actual phosphorylation sequence
# the alignment, and the parsed phosphosite array
# Output true or false, does the site given in the file match the neighboring sequence
def TestPhosphosites(neighborhood, phosphosites, Alignment):

    for n in range(0,19):#len(neighborhood)):
        for m in range(0,len(neighborhood[n])):
            pos = FindinList('p', neighborhood[n][m], m+1)
            nump = 0
            for i in range(0,pos):
                if i<pos-nump:
                    if neighborhood[n][m][i+nump] == 'p':
                        nump = nump+1
                    if Alignment[n][0][TranslatePos2Align(phosphosites[n][m]-(pos-i), Alignment[n][0])] != neighborhood[n][m][i+nump]:
                        print "FALSE" 
                        print neighborhood[n][m]
                        print neighborhood[n][m][i+nump]
                        print i
                        print ('phosphosite', phosphosites[n][m], 'phosphobeginning', phosphosites[n][m]-(pos-i), i+nump)
                        print Alignment[n][0][TranslatePos2Align(phosphosites[n][m]-(pos-i), Alignment[n][0])]
                        print nump
                        print ('M', m, HowManySites(neighborhood[n][m]))
                        print Alignment[n][0][TranslatePos2Align(phosphosites[n][m]-5-(pos-i),Alignment[n][0]):TranslatePos2Align(phosphosites[n][m]+5-(pos-i), Alignment[n][0])]
                        print '\n'
#for i in range(pos,len(neighborhood[n][m])):
#               if Alignment[n][0][TranslatePos2Align(phosphosites[n][m]+1+i, Alignment[n][0])] != neighborhood[n][m][i]:
#                   print false      


#Input: tab delimited files with phosphorylation site position in S. Cer for every gene
# filename is: './PhosphorylationSites.tab
def ParsePhosphoSiteFilesNeighborhood(filename, genes):

    SCerPhosphoNeighborhood = []
    
    for i in range(0, len(genes)):
        SCerPhosphoNeighborhood.append([])
        File = open('./' + filename, 'r')
        
    line = File.readline()
    while (line):
        line = File.readline()
        geneindex = FindinList(str(line.split('\t')[0]).upper(),genes)
            
        if geneindex != -1 and line != '':
             SCerPhosphoNeighborhood[geneindex].append(str(line.split('\t')[2]))
        else:
             print line.split('\t')[0]

    return SCerPhosphoNeighborhood


def TestGenes(ga):
    badgenes=[]
    for i in range(0, len(ga[1])):   
        count = 0
        for j in range(0, len(ga[1][i])):
            if ga[1][i][j] == []:
                count = count+1
            if count == len(ga[1][i]):
                badgenes.append(ga[0][i])
    badgenes
    for i in range(0, len(badgenes)):
        shutil.move('./Alignments/'+str(badgenes[i])+'.aln', './badgenes/')


    return badgenes

def HowManySites(seq):
    count = 0
    for i in range(0,len(seq)):
         if seq[i] == 'p':
            count = count+1
    
    return count


