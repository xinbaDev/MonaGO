import sys
sys.path.append('../')

import logging
import traceback as tb
import suds.metrics as metrics
from suds import *
from suds.client import Client
from datetime import datetime

errors = 0

# setup_logging()

# logging.getLogger('suds.client').setLevel(logging.DEBUG)

url = 'http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl'
    
print 'url=%s' % url

#
# create a service client using the wsdl.
#
client = Client(url)
client.wsdl.services[0].setlocation('https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')
#
# print the service (introspection)
#
#print client

#authenticate user email 
client.service.authenticate('zxin9@student.monash.edu')

#add a list 
inputIds = '31741_at,31734_at,32696_at,35957_at,39519_at'
idType = 'AFFYMETRIX_3PRIME_IVT_ID'
listName = 'make_up'
listType = 0
client.service.addList(inputIds, idType, listName, listType)

#print client.service.getDefaultCategoryNames()

client.service.setCategories("GOTERM_BP_FAT");


print client.service.getTableReport()
# print client.service.getGeneReportCategories();
#getTermClusterReport
# overlap=3
# initialSeed = 3
# finalSeed = 3
# linkage = 0.5
# kappa = 20 
# print client.service.getTermClusterReport(overlap, initialSeed, finalSeed, linkage, kappa)



