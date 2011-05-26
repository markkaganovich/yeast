class kwalalignindex:
    args = 'kwalalign'
    fun = 'getalignindex'

class kwalallaacons:
    args = ['kwalalign', 'seqpos']
    fun = 'seqconservation'

class kwalpsitecons:
    args = ['kwalalign', 'phosphosites']
    fun = 'seqconservation'

class kwalsytcons:
    args = ['kwalalign', 'sytpos']
    fun = 'seqconservation'

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

class syt:
    fun = 'getSYTpos'
    
class sytturnover:
    args = ['sytturnover',
                          'kwalalign',
                          'kwalalignindex',
                          'Scer',
                          'Kwal',
                          [0, 5, 10],
                          'sytpos']
    fun = 'calcturnover'


