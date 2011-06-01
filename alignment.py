import os

def align(alignmentprogram = 'mafft', inputpath = './data/Alignments/input/', outputpath = './data/Alignments/output/'):
    for file in os.listdir(inputpath):
        print alignmentprogram+' '+inputpath+file+' > '+ outputpath+file
        if alignmentprogram == 'muscle':
            os.system('~/Downloads/muscle' + ' -in ' + inputpath + file + ' -out ' + outputpath+file)
        else:
            os.system(alignmentprogram+' '+inputpath+file+' > '+ outputpath+file)

def filestodic(alignmentfiles = './data/Alignments/ScerKwal/'):
    alignment ={}
    for f in os.listdir(alignmentfiles):
        FILE = open(alignmentfiles+f)
        lines = FILE.readlines()
        f = f.split('.')[0]
        alignment[f] = {}
        for line in lines:
            if line.startswith('>'):
                species = line.strip('\n').split('>')[1]
                seq = ''
            else:
                seq = seq + line.strip('\n')
            alignment[f][species] = seq

    return alignment

def makefastainput(ref_species, orthospecies, orthologdic, **sequencefiles):
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
#os.mkdir('./data/Alignments/' + orthospecies)
#       os.mkdir('./data/Alignments/' + orthospecies + '/input/')
        file = open('./data/Alignments/' + orthospecies +'/input/'+g+'.fasta','w')
        file.write('>'+ref_species+'\n'+ref[g]+'\n')
        for i in sequencefiles.keys():
            if not i == ref_species:
                if g in orthologdic.keys():
                    orthogene = orthologdic[g]
                if not orthogene == 'NONE' and orthogene in genomes[i].keys():
                    file.write('>'+i+'\n'+genomes[i][orthogene]+'\n')
        file.close()

'''
align two sequences
'''
def alignseqs(seq1, seq2, name1, name2):
    if not 'tempseqs' in os.listdir('./'):
        os.mkdir('./tempseqs')
    filename = name1 + '_' + name2 + '.fasta'
    file = open('./tempseqs/' + filename,'w')
    file.write('>'+name1 + '\n' + seq1 + '\n')
    file.write('>'+name2 + '\n' + seq2)
    file.close()
    print name1+name2
    print ('mafft'+' '+ './tempseqs/' + filename +' > '+ './tempoutputs/'+filename)
    os.system('mafft ./tempseqs/' + filename +' > ./tempoutputs/'+ filename)



        



    









