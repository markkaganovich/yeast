from numpy import *
from helpers import *
import YeastPhosphoProject 
from scripts import *
import simplejson as json 


#alignment = json.load(open('alignmentfile'))
speciesalign = json.load(open('data/speciesalign'))
scergenes = json.load(open('data/scergenes'))
#phosphosites = json.load(open('phosphosites'))
orthophosphosites = json.load(open('data/orthophosphosites'))
#genenames = json.load(open('genenamesfile'))

#pairs that occurred in the whole genome duplication 
#source: 'proteins' folder from Kellis et al. paper
wgdpairs = json.load(open('wgdpairs'))

#paralogs... also from that paper... same common ancestor Kwal
paralogs = json.load(open('paralogs'))

#the intersection of these two:
wgdpair1 = json.load(open('wgd_paralog_pair1'))
wgdpair2 = json.load(open('wgd_paralog_pair2'))

wgd = json.load(open('wgd_paralogs'))
awgd = json.load(open('after_wgd_paralogs'))
bwgd = json.load(open('before_wgd_paralogs'))


#alignment of pairs of paralogs - generated by ParParse in YeastPhosphoConservation
paralignment = json.load(open('parsalignment'))
parphosphosites = json.load(open('parphosphosites'))



orthos = GenOrthologs(scergenes)
DICparsalignment = ParParse(scergenes, 'WGD_pars_alignments')
parphosphosites_dic= ParseParPhosphoSiteFiles('./PhosphorylationSites.tab', scergenes, DICparsalignment, orthos)











