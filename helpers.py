# Helper functions for [yeast] project

import YeastPhosphoProject 
import sys
import os, subprocess

# Wapinski supplementary file... paralog pairs and their origin
def ParseParalogFile(list):
    parWGD, parA, parB, parC, parD, parG, parE, parH, parI, parJ =[],[], [],[],[],[],[],[],[],[]
    for i in range(2,len(list)):
        a = list[i].split('\t')
        eval('par'+str(a[2])).append([a[0],a[1]])
    
    return [parWGD, parA, parB, parC, parD, parE, parG, parH, parI, parJ]

#### convert ParseParalogFile to inds
def Parinds(parGenes, scergenes):
    Parind=[]
    for par in parGenes:
        if par[0] in scergenes and par[1] in scergenes:
            Parind.append([scergenes.index(par[0]), scergenes.index(par[1])])

    return Parind




def GeneList(filename = './orf_trans.fasta'):
    genes = []
    File = open(filename)
    line = File.readline()
    while(line):
        if line[0] == '>':
            words = line.split('\n')
            genes.append(words[0][1:len(words[0])])
        line = File.readline()

    savedfile = filename.split('.')[-2].split('/')[-1]
    SavetoFile(genes, '../genomes/'+savedfile+'Genes')

    return genes    

def GenSpecies():

     species = ['Spar', 'Smik', 'Skud', 'Sbay', 'Scas', 'Sklu', 'Klac', 'Calb', 'Cgla',  'Spom']

     return species
    
# pairwise species paralog inpus
def makeparalignmentinput(orthos,pars):

    species = GenSpecies()
    
    

    for s in range(0,len(species)):
        if not species[s] == 'Skud':
            genenames = LoadfromFile('../genomes/'+species[s]+'Genes')
            geneseq = LoadfromFile('../genomes/'+species[s]+'Seq')
            for i in range(0, len(pars)):
                print species[s]
                print orthos[pars[i][0]][s]
                
                if  not orthos[pars[i][0]][s] == [] and not orthos[pars[i][1]][s] == []:
                    file = open('../genomes/'+species[s]+'_paraligninputs/'+str(genenames[orthos[pars[i][0]][s]])+'_'+str(genenames[orthos[pars[i][1]][s]]), 'w')
                    file.write('>'+str(genenames[orthos[pars[i][0]][s]]) +'\n')
                    file.write(str(geneseq[orthos[pars[i][0]][s]])+'\n')
                    file.write('>'+str(genenames[orthos[pars[i][1]][s]]) +'\n')
                    file.write(str(geneseq[orthos[pars[i][1]][s]])+'\n')
                


def GenOrthologs(genes):
    
    species = ['Spar', 'Smik', 'Skud', 'Sbay', 'Scas', 'Sklu', 'Klac', 'Calb', 'Cgla',  'Spom']

    Orthologs = []
    for o in range(0,len(genes)):
        Orthologs.append([])
        for s in range(0,len(species)):
            Orthologs[o].append([])
    for i in range(0,len(species)):
        print i
        if i==2:
            i=i+1
            continue 
        File = open('./data/orthologfiles/Scer-'+str(species[i])+'-orthologs.txt')
        genes2 = LoadfromFile('../lab/genomes/'+species[i]+'Genes')

        lines = File.readlines()
        for line in lines:
            gene = str(line).split('\t')
            index = FindinList(str(gene[0]), genes)
            par = str(gene[1])
            if i ==0 or i==1:
                par = par[0:4]+'_'+par[4:len(par)]
            if i==9:
                par = par[0:4].upper() + par[4:len(par)]
            if i==4:
                temp = par.split('\n')
                par = temp[0]
            index2 = FindinList(par, genes2)
            if index== -1:
                print gene[0]
            else:
                Orthologs[index][i] = index2

    return Orthologs

# Sept 10, 2009 Find Paralogs of each Scer gene and their respective orthologs in the other species    
def GenParalogs():

    genes= LoadfromFile('../genomes/ScerGeneNames')

    Paralogs = []
    for i in range(0, len(genes)):
        Paralogs.append([])

    p = LoadfromFile('../genomes/Scer-paralogs.txt')
    for i in range(0,len(p)):
        if p[i].split('\t')[1] != 'NONE':
            g = FindinList(p[i].split('\t')[0], genes) 
            if g!=-1:
                Paralogs[g] = FindinList(p[i].split('\t')[1], genes)

    return Paralogs
            
# process Paralogs

# Check for ones that have the same ortholog in Klac ... this narrows it down to a group of genes that likely had a definite common ancestral gene that was duplicated
def OneOrtho(paralogs, orthos):

    klac = 6
    singlepars = [-1]*len(paralogs)
    pargenes = []

    for i in range(0,len(paralogs)):
        if paralogs[i] != []:
            if orthos[i][klac] == orthos[paralogs[i]][klac]:
                singlepars[i] = paralogs[i]
            singlepars[paralogs[i]] = -1

    SavetoFile(singlepars, '../genomes/klacsinglepars')
    
    for i in range(0, len(singlepars)):
        if singlepars[i]!=-1:
            pargenes.extend([i])
            pargenes.extend([singlepars[i]])

    return [singlepars, pargenes]


# compare paralogs in SCer to each other
def compareparalogs(pars, genes):
 
    files = os.listdir('../genomes/AllGeneSeqs')
    pairs = []

    print files
    for i in range(0,len(pars)):
        print i
        if genes[pars[i][0]]+'.fasta' in files and genes[pars[i][1]]+'.fasta' in files:
            print genes[pars[i][0]]
            File = open('../genomes/AllGeneSeqs/'+genes[pars[i][0]]+'.fasta')
            scerseq1 = File.readline()
            File = open('../genomes/AllGeneSeqs/'+genes[pars[i][1]]+'.fasta')
            scerseq2 = File.readline()

            File = open('../genomes/WGD_parsaligninputs/'+genes[pars[i][0]]+'_'+genes[pars[i][1]]+'.fasta', 'w') 
            File.write('>' + genes[pars[i][0]] + '\n')
            File.write(scerseq1)
            File.write('>' + genes[pars[i][1]] + '\n')
            File.write(scerseq2)
            File.close

#  print "../genomes/parsaligninputs/" + genes[i] +'-'+genes[pars[i]] 

#os.system("mafft "+ "../genomes/parsaligninputs/YAL005C-YER103W"+ ">"+ 'testalignYAL005C' )
#           result= subprocess.call(["mafft", '../genomes/parsaligninputs/YAL005C-YER103W'])
#           print result
#           result = subprocess.call(["mafft", '../genomes/parsaligninputs/'+genes[i]+'-'+genes[pars[i]]])
#           print result

#os.system("mafft "+ "../genomes/parsaligninputs/" + genes[i] +'-'+genes[pars[i]]+" >"+" ../genomes/parsalignoutput/"+genes[i]+'-'+genes[pars[i]]+".output")



    

# generate inputs for MAFFT, fasta format with one file per gene and each file
# contains the orthologs of the gene in the other species
def GenInputs(orthos, ScerGenes):

    species = GenSpecies()
    Seq = []

    for s in range(0, len(species)):
        Seq.append(LoadfromFile('../genomes/'+species[s]+'Seq'))
        
    
    for i in range(0, len(ScerGenes)):
        if orthos[i][0]!=[] and orthos[i][0]!=-1:
            File = open('../genomes/AlignmentInput2/'+str(ScerGenes[i])+'.fasta', 'w')
            File.write('>'+'Scer')
            File.write('\n')
            File.write(Seq[0][i])
            File.write('\n')
            for j in range(1, len(species)):
                if orthos[i][j-1] != [] and orthos[i][j-1] != -1:
                    File.write('>'+species[j])
                    File.write('\n')
                    File.write(Seq[j][orthos[i][j-1]])
                    File.write('\n')

#parse orf-trans-all.fasta file
def ParseOrfs(orffile, scergenes):
    file = open(orffile)
    lines = file.readlines()
    seq=''
    for i in range(0, len(lines)):
        if lines[i][0] == '>':
            word =lines[i].split(' ')[0]
            if seq != '':
                out = open('../genomes/AllGeneSeqs/'+genename+'.fasta', 'w')
                out.write(seq + '\n')
                out.close()
            seq=''
            genename = word[1:len(word)]
        else:
            if genename in scergenes:
                seq = seq+lines[i][0:-1]

def get_gene_synonyms(orffile, scergenes):
    file = open(orffile)
    lines = file.readlines()
    syn = []
    for i in range(0, len(scergenes)):
        syn.append([])

    for i in range(0, len(lines)):
        if lines[i][0] == '>':
            genename = lines[i].split(' ')[0]
            genename = genename[1:len(genename)]
            synword = lines[i].split(' ')[1] 
            if genename in scergenes:
                syn[scergenes.index(genename)] = synword


    return syn

    
def GenSequences():
    
    species = GenSpecies()

    for i in range(0, len(species)):
        File = open('../genomes/'+str(species[i])+'.fasta')

        genes = LoadfromFile('../genomes/'+str(species[i])+'Genes')
        
        seq = []
        for j in range(0,len(genes)):
            seq.append('')


        line = 1
        f=0
        while(line):
            f = f+1
            line = File.readline()
            name = line.split('\n')[0]
            name = name[1:len(name)]
            index = FindinList(name, genes)
            line = File.readline()
            line = line.split('\n')[0]
            if index != -1 and name!='':
                seq[index] = line
       
        SavetoFile(seq, '../genomes/'+species[i]+'Seq')







def SavetoFile(variable, filename):
    File = open(filename, 'w')

    for i in range(0,len(variable)):
        if isinstance(variable[i], type(0)):
            File.write(str(variable[i]))
            File.write('\n')
        else:    
            File.write(variable[i])
            File.write('\n')

def LoadfromFile(filename):
    File = open(filename)

    variable =[]
    line = File.readline()
    while(line):
        variable.append(line.split('\n')[0])
        line = File.readline()

    return variable

def getridofnan(list):
    return filter(lambda x: str(x) != 'nan', list)
        
##############################################################
##############################################################
##############################################################
##############################################################

# basically just need to run this once
def ParseFastaGenes():
    SparGenes = GeneList('../genomes/Spar.fasta')
    SmikGenes = GeneList('../genomes/Smik.fasta')
    SbayGenes = GeneList('../genomes/Sbay.fasta')
    ScasGenes = GeneList('../genomes/Scas.fasta')
    KlacGenes = GeneList('../genomes/Klac.fasta')
    SkluGenes = GeneList('../genomes/Sklu.fasta')
    KwalGenes = GeneList('../genomes/Kwal.fasta')
    CalbGenes = GeneList('../genomes/Calb.fasta')
    SpomGenes = GeneList('../genomes/Spom.fasta')

# generate gene names from ScerGenes files (the ScerGenes, etc. files have names and descriptions)
def GenGeneNames():
    genes =[]
    species = GenSpecies()
    for i in range(0, len(species)):
        fullgenes = LoadfromFile('../genomes/'+species[i]+'Genes')
        for j in range(0, len(fullgenes)):
            genes.append(fullgenes[j].split('t')[0])
        SavetoFile(genes, '../genomes/'+species[i]+'GeneNames')

## how many genes does each sequences have
def howmanygenes():
    species = GenSpecies()
    for i in range(0, len(species)):
        count = 0
        File= open('../genomes/'+species[i]+'Seq')

        line  = 1
        while(line):
            line = File.readline()
            if len(line)>10:
                count = count+1
            
        print species[i]
        print count
        




