import simplejson

datasource = '/Users/markkaganovich/Dropbox/data/'

def json(filename, datasource =''):
    file = open(datasource + filename)
    result = simplejson.load(file)
    file.close()
    return result


