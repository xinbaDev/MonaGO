from HTMLParser import HTMLParser
import re
import pycurl
from StringIO import StringIO
import certifi
import logging
import time
import select

from twisted.internet import defer

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

'''
pycurl helper class

'''
class PycurlHelper:

    def __init__(self):
        self.curl = pycurl.Curl()
        self.curls = []#for multicurl
        self.buffers = []
        self.i = 0

        self.s = pycurl.CurlShare()
        self.s.setopt(pycurl.SH_SHARE, pycurl.LOCK_DATA_COOKIE)
        self.s.setopt(pycurl.SH_SHARE, pycurl.LOCK_DATA_DNS)
        self.s.setopt(pycurl.SH_SHARE, pycurl.LOCK_DATA_SSL_SESSION)

    def get(self,url):
        buffer = StringIO()
        self.curl.setopt(pycurl.URL, url)
        self.curl.setopt(pycurl.CUSTOMREQUEST, "GET")
        self.curl.setopt(self.curl.WRITEDATA, buffer)
        self.curl.perform()
        return buffer.getvalue().decode('iso-8859-1')

    def post(self,url,data):
        return 0


    def sendMultipart(self,url,data):
        buffer = StringIO()
        self.curl.setopt(pycurl.URL, url)
        self.curl.setopt(pycurl.HTTPPOST, data)
        self.curl.setopt(pycurl.SHARE, self.s)
        self.curl.setopt(pycurl.CUSTOMREQUEST, "PUT")
        self.curl.setopt(self.curl.WRITEDATA, buffer)
        self.curl.setopt(pycurl.COOKIEJAR, 'cookie.txt')
        self.curl.setopt(pycurl.CAINFO, certifi.where())

        start_time = time.time()
        self.curl.perform()
        logger.info("upload data lasts--- %s seconds ---" % (time.time() - start_time))
        
        return buffer.getvalue().decode('iso-8859-1')

    def close(self):
        self.curl.close()


    def CurlMultiGet(self,urls):

        SELECT_TIMEOUT = 5.0

        m = pycurl.CurlMulti()

        curls = []

        for url in urls:
            curl = pycurl.Curl()
            #curl.setopt(pycurl.VERBOSE, 1)#for debugging, show details about connection status 
            curl.setopt(pycurl.SHARE,self.s)#share ssl session,cookies,DNS cache
            curl.setopt(pycurl.URL, url)
            curl.setopt(pycurl.CUSTOMREQUEST, "GET")

            mybuffer = StringIO() 
            curl.setopt(pycurl.WRITEDATA, mybuffer)
            curl.setopt(pycurl.CAINFO, certifi.where())

            self.buffers.append(mybuffer)
            curls.append(curl)

        #m.add_handle does not add a Python reference to the Curl object
        #therefore need a container to hold the reference to curl object
        for curl in curls:
            m.add_handle(curl)         

        start_time = time.time()

        while True:
            ret, num_handles = m.perform()
            if  ret != pycurl.E_CALL_MULTI_PERFORM:
                break

        while num_handles:

            ret = m.select(SELECT_TIMEOUT)

            if ret == -1:
                continue
            while 1:
                ret, num_handles = m.perform()

                if ret != pycurl.E_CALL_MULTI_PERFORM: 
                    break

        
        num_q, ok_list, err_list = m.info_read()

        
        logger.info("run pycurl multi lasts %s" % (time.time() - start_time))

        # for b in self.buffers:
        #         logger.debug("buffer:{}".format(b.getvalue()))




class DavidDataScrawler(object):

    def setParams(self,inputIds,idType,annotCat,pVal):
        self.inputIds,self.idType,self.annotCat,self.pVal = inputIds,idType,annotCat,pVal

    def run(self):

        pcHelper = PycurlHelper()

        #d = Deferred()#init defer

        res = self._uploadGene(pcHelper,self.idType,self.inputIds)

        if self._checkSuccess(res):

            url_1 = 'https://david.ncifcrf.gov/chartReport.jsp?annot={0}&currentList=0'.format(self.annotCat)
            url_2 = 'https://david.ncifcrf.gov/list.jsp'

            urls = [url_1,url_2]

            pcHelper.CurlMultiGet(urls)

            getGO_response,geneList_response = map(lambda x: x.getvalue().decode('iso-8859-1'), pcHelper.buffers)
            
            #logger.debug("getGO_response:{}".format(getGO_response))
            #logger.debug("geneList_response:{}".format(geneList_response))

            go,geneIds = self._parseGO(getGO_response, pcHelper)
            geneList = self._parseGenes(geneList_response)
                

            go_filtered = self._filterGO(self.pVal,go)
            geneIds = self._getUniqueGeneIds(geneIds)

            #logger.debug("geneIds:{}".format(geneIds))

            geneIdNameMapping = self._getGenesNamesByIds(geneIds,geneList)

            #change the gene id into gene name in go
            go_processed = self._changeGeneIdToNameInGO(go_filtered,geneIdNameMapping)

            if not go_processed:
                raise Exception("get final GO failed")
                

            return go_processed


        else:
            logger.info("get chartReport failed")
            raise Exception("upload genes to DAVID failed")


    def _checkSuccess(self,res):
        if res.find("DAVID: Functional Annotation Tools")==-1:
            return False
        else:
            return True


    def _parseGO(self,getGO_response,pcHelper):
        parser = DavidDataScrawler.GOParser(pcHelper)
        parser.feed(getGO_response)#get go
        go = parser.getGO_inf()
        geneIds = parser.getGeneLists()

        #logger.debug('go:{}'.format(go))

        if not go:
            raise Exception("get GO terms failed") 

        #logger.debug('geneIds:{}'.format(geneIds))

        if not geneIds:
            raise Exception("get gene lists failed") 

        return go,geneIds


    def _parseGenes(self,geneList_response):
        parser = DavidDataScrawler.geneParser()
        parser.feed(geneList_response)
        genesIdName = parser.getParsedGeneNameWithId()
        return genesIdName





    def _uploadGene(self,pcHelper,idType,inputIds):

        data = [('idType', idType), ('uploadType', 'list'),('multiList','false'),('Mode','paste'),
                         ('useIndex','null'),('usePopIndex','null'),('demoIndex','null'),('ids',inputIds),
                         ('removeIndex','null'),('renameIndex','null'),('renamePopIndex','null'),('newName','null'),
                         ('combineIndex','null'),('selectedSpecies','null'),('uploadHTML','null'),('managerHTML','null'),
                         ('sublist',''),('rowids',''),('convertedListName','null'),('convertedPopName','null'),
                         ('pasteBox',inputIds),('Identifier',idType) , ('rbUploadType','list')]

        return pcHelper.sendMultipart(url="https://david.ncifcrf.gov/tools.jsp",data=data)




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
            if len(geneIds) > 1:
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

        def __init__(self,pcHelper):
            HTMLParser.__init__(self)

            self.pcHelper = pcHelper
            self.go_inf = []
            self.geneLists = {}
            self.metacount = 0

        def handle_starttag(self, tag, attrs):

            #get GO_id,GO_name,p-value,count
            if tag == "a":
                m = re.search('(data/download/chart_\w+.txt)',attrs[0][1])
                if m!=None:
                    url = 'https://david.ncifcrf.gov/'+m.group(0)
                    res = self.pcHelper.get(url)
                    self._parseGO(res)

            #get gene rowid
            if tag == "img":
                if attrs[0][1] == 'graphics/two_tone_2_a.jpg':
                    genes = attrs[6][1].split(";")[1]
                    self.geneLists[self.metacount] = self._parseStringIntoList(genes)
                    self.metacount+=1

        def _parseGO(self,GO):
            line = GO.encode('ascii','ignore').split("\n")
            for i in range(1,len(line)-1):
                cell = line[i].split("\t")
                GO_term = cell[1].split("~")
                self.go_inf.append({"cat":cell[0],"GO_id":GO_term[0],"GO_name":GO_term[1],"count":cell[2],"pVal":cell[4],"genes":cell[5].lower().strip()})

        def getGO_inf(self):
            return self.go_inf

        def getGeneLists(self):
            return self.geneLists

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