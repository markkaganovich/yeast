import os

def align(alignmentprogram = 'mafft', inputpath = './data/Alignments/input/',
        outputpath = './data/Alignments/output/'):
    for file in os.listdir(inputpath):
        print alignmentprogram+' '+inputpath+file+' > '+ outputpath+file
        os.system(alignmentprogram+' '+inputpath+file+' > '+ outputpath+file)

def filestodic(alignmentfiles = './data/Alignments/ScerKwal/'):
    alignment ={}
    for file in os.listdir(alignmentfiles):
        FILE = open(alignmentfiles+file)
        lines = FILE.readlines()
        if lines:
            genename = file.split('.')[0]
            alignment[genename] = {}
            for line in lines:
                if line.startswith('>'):
                    species = line.split('>')[1]
                    species = species.strip('\n')
                else:
                    alignment[genename][species] = line.strip('\n')

    return alignment

def makefastainput(ref_species, orthologdic, **sequencefiles):
    
    genomes = {}
    for k in sequencefiles.keys():
        print k
        file = open(sequencefiles[k])
        lines = file.readlines()
        file.close()
        genomes[k] = {}
        for line in lines:
            if line.startswith('>'):
                orfname = line.split('>')[1]
            else:
                genomes[k][orfname.strip('\n')] = line.strip('\n')

    ref = genomes[ref_species]
    for g in ref.keys():
        file = open('./data/Alignments/input/'+g+'.fasta','w')
        file.write('>'+ref_species+'\n'+ref[g]+'\n')
        for i in sequencefiles.keys():
            if not i == ref_species:
                orthogene = orthologdic[i][g]
                if not orthogene == 'NONE':
                    file.write('>'+i+'\n'+genomes[i][orthogene]+'\n')
        file.close()




        



    









