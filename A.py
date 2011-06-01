class seq:
    fun = 'getseq'
    data = 'ScerSeqs'
    
class init:
    fun = '__init__'
    data = 'scergenes'

class phosphositeturnover:
    args = ['phosphositeturnover', 'kwalalign', 'kwalalignindex', 'Scer',
'Kwal', [0, 5, 10], 'phosphosite']
    fun = 'calcturnover'

class scasalignindex:
    args = 'scasalign'
    fun = 'getalignindex'

class scasallaacons:
    args = ['scasalign', 'seq']
    fun = 'seqconservation'

class spomalignindex:
    args = 'spomalign'
    fun  = 'getalignindex'

class sytturnover:
    args = ['sytturnover',
                          'kwalalign',
                          'kwalalignindex',
                          'Scer',
                          'Kwal',
                          [0, 5, 10],
                          'sytpos']
    fun = 'calcturnover'

class calbalign:
    fun = 'getotheralign'
    args = 'calbalign'
    data = 'Calbalign'

class orthologs:
    fun = 'getorthos'
    data = 'KellisOrthos'
    
class kwalalign:
    fun = 'getalign'
    data = 'ScerKwalAlignment'

class goterms:
    fun = 'getGOterms'
    data = 'goterms'

class speciesalign:
    fun = 'getmultiplealign'
    data = 'speciesalign'

class genename:
    fun = 'genename'
    data = 'genenames'

class scasalign:
    fun = 'getotheralign'
    args = 'scasalign'
    data = 'scasalign'

class spomalign:
    fun = 'getotheralign'
    args = 'spomalign'
    data = 'spomalign'

class histonemod:
    fun = 'hasproperty'
    args = 'histonemod'
    data = 'histonemodifiersJSON'

class TF:
    fun = 'hasproperty'
    args = 'TF'
    data = 'TFsJSON'

class phosphosites:
    fun = 'getphosphosites'
    data = 'gersteinphosphositefile'

class wapinskiparalogs:
    fun = 'getWparalogs'
    data = 'paralogslist'

class initparsW:
    fun = 'initall'
    data = 'plist'

class paralign:
    fun = 'getalign'
    data = 'Wparsalignment'

class alignmentstats:
    fun = 'alignmentstats'

class phosphostats:
    fun = 'phosphostats'

class event:
    fun = 'getevents'
    data = 'wapinskieventsdic'


