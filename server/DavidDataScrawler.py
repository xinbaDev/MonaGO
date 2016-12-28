from HTMLParser import HTMLParser
import re
from PycurlHelper import PycurlHelper

from cookielib import CookieJar
from twisted.internet import defer, reactor
from twisted.web.client import Agent, readBody, HTTPConnectionPool, CookieAgent 
from twisted.web.http_headers import Headers
from twisted.internet.ssl import ClientContextFactory

from uuid import uuid4
from MultiPart import MultiPartProducer

import logging
from logTime import logTime

logging.basicConfig(level=logging.debug)
logger = logging.getLogger(__name__)

class WebClientContextFactory(ClientContextFactory):
    def getContext(self, hostname, port):
        return ClientContextFactory.getContext(self)

class DavidDataScrawler(object):

    def run(self):
        

        self.go_inf = []
        self.go_processed = []

        pool = HTTPConnectionPool(reactor)
        contextFactory = WebClientContextFactory()

        cookieJar = CookieJar()
        agent = CookieAgent(Agent(reactor,contextFactory, pool=pool),cookieJar)

        d = self._uploadGene(agent,self.idType,self.inputIds)

        d.addCallback(lambda ign: self.handleResponse(agent))

        reactor.run(installSignalHandlers=0)


    @logTime
    def _uploadGene(self,agent,idType,inputIds):

        data = [('idType', idType), ('uploadType', 'list'),('multiList','false'),('Mode','paste'),
                         ('useIndex','null'),('usePopIndex','null'),('demoIndex','null'),('ids',inputIds),
                         ('removeIndex','null'),('renameIndex','null'),('renamePopIndex','null'),('newName','null'),
                         ('combineIndex','null'),('selectedSpecies','null'),('uploadHTML','null'),('managerHTML','null'),
                         ('sublist',''),('rowids',''),('convertedListName','null'),('convertedPopName','null'),
                         ('pasteBox',inputIds),('Identifier',idType) , ('rbUploadType','list')]

        boundary = uuid4().hex
        boundaryHeader = "multipart/form-data; boundary="+boundary

        postBody = MultiPartProducer(data, boundary = boundary)
        

        d = agent.request(
            'POST',
            'https://david.ncifcrf.gov/tools.jsp',
            Headers({      
                "Content-Type": [boundaryHeader],
                }),
            postBody)

        d.addCallback(self.cbUploadGeneRequest)

        return d

    def setParams(self,inputIds,idType,annotCat,pVal):
        self.inputIds,self.idType,self.annotCat,self.pVal = inputIds,idType,annotCat,pVal

    def cbGetGeneList(self,response):
        d = readBody(response)
        d.addCallback(self.cbParseGeneList)
        return d

    def cbParseGeneList(self,geneList_response):
        print "get gene list done"
        self.geneList = self._parseGenes(geneList_response)

    def cbResponseAfterUploadGenes(self,body):
        self.res = body

    def cbUploadGeneRequest(self,response):
        d = readBody(response)
        d.addCallback(self.cbResponseAfterUploadGenes)
        return d


    def cbParseGO(self,getGO_response):
        self.url, self.geneIds = self._parseGO(getGO_response)
        print "get GO url,geneid done"


    def cbGetGOChartHtml(self,response):
        d = readBody(response)
        d.addCallback(self.cbParseGO)
        return d

    def getDeferedGOChartResponse(self,agent,url):
        d = agent.request(
            'GET', url,
            Headers({'User-Agent': ['Chrome']}),
            None)
        d.addCallback(self.cbGetGOChartHtml)
        return d


    def getDeferedGOMappingResponse(self,agent,url):
        
        d = agent.request(
            'GET', url,
            Headers({'User-Agent': ['Chrome']}),
            None)
        d.addCallback(self.cbGetGeneList)

        return d


    def cbGenerateGOData(self,res):
        print "process data"
        go_filtered = self._filterGO(self.pVal,self.go_inf)
        
        # with open("go_filtered","w") as fw:
        #     fw.write(str(go_filtered))
        geneIds = self._getUniqueGeneIds(self.geneIds)

        #logger.debug("geneIds:{}".format(geneIds))

        geneIdNameMapping = self._getGenesNamesByIds(geneIds,self.geneList)

        #change the gene id into gene name in go
        self.go_processed = self._changeGeneIdToNameInGO(go_filtered,geneIdNameMapping)


    @logTime
    def cbParseGOtxt(self,GO):
        line = GO.encode('ascii','ignore').split("\n")
        for i in range(1,len(line)-1):
            cell = line[i].split("\t")
            GO_term = cell[1].split("~")
            self.go_inf.append({"cat":cell[0],"GO_id":GO_term[0],"GO_name":GO_term[1],"count":cell[2],"pVal":cell[4],"genes":cell[5].lower().strip()})
        #might need to callback


    def cbGetGOtxt(self,response):
        d = readBody(response)
        d.addCallback(self.cbParseGOtxt)
        return d


    def cbGetGenes(self,agent):
        d = agent.request(
            'GET', self.url,
            Headers({'User-Agent': ['Chrome']}),
            None)
        d.addCallback(self.cbGetGOtxt)
        return d

    def handleResponse(self,agent):
        if self._checkSuccess(self.res):
            print "upload sucess"

            url_1 = 'https://david.ncifcrf.gov/chartReport.jsp?annot={0}&currentList=0'.format(self.annotCat)
            url_2 = 'https://david.ncifcrf.gov/list.jsp'
            deferList = []
            d1 = self.getDeferedGOChartResponse(agent,url_1)
            d1.addCallback(lambda ign: self.cbGetGenes(agent))
            d2 = self.getDeferedGOMappingResponse(agent,url_2)
            deferList.append(d1)
            deferList.append(d2)
            d = defer.DeferredList(deferList)
                
            d.addCallback(self.cbGenerateGOData)

        else:
            logger.info("get chartReport failed")
            raise Exception("upload genes to DAVID failed")


    def _checkSuccess(self,res):
        if res.find("DAVID: Functional Annotation Tools")==-1:
            return False
        else:
            return True


    def _parseGO(self,getGO_response):
        parser = DavidDataScrawler.GOParser()
        parser.feed(getGO_response)#get go
        geneIds = parser.getGeneLists()
        GOUrl = parser.getGOUrl()
        #logger.debug('geneIds:{}'.format(geneIds))

        if not GOUrl:
            raise Exception("get GOUrl failed") 

        if not geneIds:
            raise Exception("get gene lists failed") 

        return GOUrl,geneIds


    def _parseGenes(self,geneList_response):
        parser = DavidDataScrawler.geneParser()
        parser.feed(geneList_response)
        genesIdName = parser.getParsedGeneNameWithId()
        return genesIdName



    def _changeGeneIdToNameInGO(self,go,geneIdNameMapping):
        '''
        transfrom the gene ids in GO terms to gene name according to the mapping

        Args:
            go:a list of GO terms
            geneIdNameMapping:mapping between gene id and gene name

        return:
            processed GO terms
        '''

        for i in range(0,len(go)):

            geneNames = go[i]["genes"].split(",")
            geneNameString = "|".join([geneIdNameMapping[geneName.strip().lower()] for geneName in geneNames]) 

            go[i]["genes"] = geneNameString[:-1]# strip the last '|'

        return go


    def _getGenesNamesByIds(self,geneIds,genesIdName):
        '''
        create a dict(geneIdNameMapping) to map gene id to gene name

        Args:
            pcHelper: pycurl helper object
            geneIds: an array of gene ids

        return:
            a dict mapping gene id to gene name
        '''

        geneIdNameMapping = {}

        for i in range(0,len(genesIdName)-1,2):#ugly way, need improve
            geneIds = genesIdName[i].split(",")
            if len(geneIds) > 1:#if multiple gene ids
                for j in geneIds:
                    geneIdNameMapping[j.strip().lower()] = genesIdName[i+1]
            else:
                geneIdNameMapping[genesIdName[i].lower()] = genesIdName[i+1]

        return geneIdNameMapping

    def _getUniqueGeneIds(self,geneListsDic):
        rowids = set([])
        rowidstr = ""

        #map(rowids.add, [genes for genes in geneListsDic.values()]) #list is not hashable
        for genes in geneListsDic.values():
            if isinstance(genes,list):
                map(rowids.add,[gene for gene in genes])
            else:
                rowids.add(genes)

        rowidstr = ','.join(rowids)

        return rowidstr



    def _filterGO(self,pVal,go_inf):

        filterGO_inf = [GO for GO in go_inf if float(GO["pVal"]) < float(pVal)]

        logger.debug('go_filtered:{}'.format(filterGO_inf))

        if not filterGO_inf:
            raise Exception("get go_filtered failed")

        return filterGO_inf



    class GOParser(HTMLParser):
        '''
        parse a list of GO terms(name and a list of associated gene ids)
        '''

        def __init__(self):
            HTMLParser.__init__(self)
            self.go_inf = []
            self.geneLists = {}
            self.metacount = 0

        def handle_starttag(self, tag, attrs):

            #get GO_id,GO_name,p-value,count
            if tag == "a":
                m = re.search('(data/download/chart_\w+.txt)',attrs[0][1])
                
                if m!=None:
                #if exists
                    url = 'https://david.ncifcrf.gov/'+m.group(0)
                    self.url = url

            #get gene rowid
            if tag == "img":
                if attrs[0][1] == 'graphics/two_tone_2_a.jpg':
                    genes = attrs[6][1].split(";")[1]
                    self.geneLists[self.metacount] = self._parseStringIntoList(genes)
                    self.metacount+=1

        def getGeneLists(self):
            return self.geneLists

        def getGOUrl(self):
            return self.url

        def _parseStringIntoList(self,genes):
            index = genes.find("geneReport")
            geneList = genes[index+10:][2:-2].encode('ascii','ignore').split(",")
            return geneList


    class geneParser(HTMLParser):
        '''
        parse gene id and gene name
        '''

        table = 0
        tr = 0
        td = 0
        genesIdName = []
        # function to handle an opening tag in the doc

        def handle_starttag(self, tag, attrs):
            if tag == "table":
                if len(attrs) > 0:
                    if attrs[0][1] == "dataTable":
                        self.table = 1

            if self.table == 1 and tag == "tr":
                if len(attrs) > 0:
                    if attrs[0][0] == "class":
                        self.tr = 1
            if self.table == 1 and self.tr == 1 and tag == "td":
                self.td+=1

        def handle_endtag(self, tag):
            if tag == "table":
                self.table = 0

            if tag == "tr" and self.table == 1:
                self.tr = 0

            if self.table == 1 and self.tr == 0:
                self.td = 0


        def handle_data(self, data):

            if self.tr == 1:
                if data != "RG" and self.td!= 4 and "\n" not in data:
                    self.genesIdName.append(data.encode('ascii','ignore'))

        def getParsedGeneNameWithId(self):
            return self.genesIdName



    # class sessionIdParser(HTMLParser):
    #     def __init__(self):
    #         HTMLParser.__init__(self)
    #         self.sessionId = 0

    #     def handle_starttag(self, tag, attrs):
    #         if tag == 'input':
    #             if attrs[1][0] == 'name' and attrs[1][1] == 'SESSIONID':
    #                 self.sessionId = attrs[2][1]

    #     def returnSessionId(self):
    #         return self.sessionId