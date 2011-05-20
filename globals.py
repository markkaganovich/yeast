import simplejson

datasource = '/Users/markkaganovich/Dropbox/data/'

def json(filename, datasource =''):
    file = open(datasource + filename)
    result = simplejson.load(file)
    file.close()
    return result

def dump(data, filename, datasource = ''):
	file = open(datasource + filename,'w')
	simplejson.dump(data, file)
	file.close()

def addfeatures(classname, feature, hash, fun, args = None, data = None):
    methods = getattr(classname, hash)
    methods[feature] = {}
    methods[feature]['fun'] = fun
    if args != None:
        methods[feature]['args'] = args
    if data != None:
        methods[feature]['data'] = data
    dump(methods, str(hash), './dbase/')
    '''
    if 'args' in kwargs.keys():
        methods['args'] = kwargs['args']
    if 'data' in kwargs.keys():
        methods['data'] = kwargs['data']
    '''



