
import sys
sys.path.append('../')

from datetime import datetime

#for praser
from HTMLParser import HTMLParser
import re

from flask import Flask,render_template,request,send_from_directory


#for pre processing the data
import DataProcessing

import pycurl
from StringIO import StringIO
import certifi

import logging


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


remote_server = False;

if(remote_server):
    root_dir = "/root/alex/myproject/"
else:
    root_dir = ""

app = Flask(__name__)
@app.route('/index', methods=['POST','GET'])
def index():

    annotCatDict = {
        'GOTERM_BP_FAT':'25',
        'GOTERM_CC_FAT':'32',
        'GOTERM_MF_FAT':'39'
    }

    if request.method == 'GET':
        return render_template("index.html")
    else:

        inputIds = request.form['inputIds']
        idType = request.form['idType']
        listName = "test"
        listType = request.form['listType']
        annotCat = request.form['annotCat']
        pVal = request.form['pVal'];

        annotCat = ','.join([annotCatDict[cat] for cat in annotCat.split(",")])

        data = getDataFromDavid(inputIds,idType,annotCat,pVal)

        matrix_count,array_order,go_hier,go_inf_reord,clusterHierData = processedData(data)

        with open(root_dir+'templates/chord_layout.html',"r") as fr_html:
            html = "".join(fr_html.readlines())

        data = "<script>"+"var go_inf="+str(go_inf_reord)+";"+"var matrix="+matrix_count+";"+"var array_order="+array_order+";"\
        +"var clusterHierData="+clusterHierData +";"+"var size="+str(len(go_inf_reord))+";"+"var goNodes="+str(go_hier)+"</script>"


    return data+html
        
        

@app.route('/css/<fileName>')
def getCss(fileName):
    return send_from_directory(root_dir+'css', fileName)

@app.route('/img/<fileName>')
def getImg(fileName):
    return send_from_directory(root_dir+'img',fileName)

@app.route('/js/<fileName>')
def getJs(fileName):
    return send_from_directory(root_dir+'js', fileName)

@app.route('/fonts/<fileName>')
def getFonts(fileName):
    return send_from_directory(root_dir+'fonts', fileName)

@app.route('/txt/<fileName>')
def getText(fileName):
    return send_from_directory(root_dir+'text', fileName)

@app.route('/demo')
def returnDemo():
    with open("visitors.txt",'a') as fw:
        fw.write("remote address: {}  time: {}\n".format(request.remote_addr,datetime.today()))
    with open(root_dir+"templates/chord_layout.html","r") as fr_html:
        html = "".join(fr_html.readlines())
    with open(root_dir+"demo/Data.txt","r") as fr:
        data = fr.readline()

    return data + html;

@app.route('/my/img/my_logo.jpg')
def getMyLogo():
    with open("visitors.txt",'a') as fw:
        fw.write("remote address: {}  time: {}\n".format(request.remote_addr,datetime.today()))
    return send_from_directory(root_dir+'my/img','my_logo.jpg')


def getDataFromDavid(inputIds,idType,annotCat,pVal):
    davidScrawler = DavidDataScrawler()
    davidScrawler.setParams(inputIds,idType,annotCat,pVal);

    try:
        go_inf_filtered_geneName = davidScrawler.run()
    except Exception as e:
        logger.error(str(e))
    else:
        return go_inf_filtered_geneName



def processedData(go_inf_filtered_geneName):

    preProcessedData = DataProcessing.dataPreprocess(go_inf_filtered_geneName)

    matrix_count = str(preProcessedData["matrix"]["matrix_count"])

    array_order = str(preProcessedData["gen_anno_reord"])

    clusterHierData = str(preProcessedData["clusterHierData"])

    go_inf_reord = preProcessedData["go_inf"]

    go_hier = preProcessedData["go_hier"]

    return matrix_count,array_order,go_hier,go_inf_reord,clusterHierData


class DavidDataScrawler(object):


    def setParams(self,inputIds,idType,annotCat,pVal):
        self.inputIds,self.idType,self.annotCat,self.pVal = inputIds,idType,annotCat,pVal

    def run(self):

        pcHelper = DavidDataScrawler.PycurlHelper()

        res = self._uploadGene(pcHelper,self.idType,self.inputIds)

        if self._checkSuccess(res):
            url = 'https://david.ncifcrf.gov/chartReport.jsp?annot={0}&currentList=0'.format(self.annotCat)
            getGO_response = pcHelper.get(url)

            parser = DavidDataScrawler.GOParser(pcHelper)
            parser.feed(getGO_response)#get go_inf
            go_inf = parser.getGO_inf()

            logger.debug('go_inf:{}'.format(go_inf))

            geneLists = parser.getGeneLists()#get a list of genes for each GO term

            logger.debug('geneLists:{}'.format(geneLists))

            go_inf_filtered = self._filterGO(self.pVal,go_inf)

            logger.debug('go_inf_filtered:{}'.format(go_inf_filtered))

            geneIds = self._getUniqueGeneIds(geneLists)

            logger.debug("geneIds:{}".format(geneIds))

            geneIdNameMapping = self._getGenesNamesByIds(pcHelper,geneIds)

            go_inf_filtered_geneName = self._changeGeneIdToNameInGO(go_inf_filtered,geneIdNameMapping)#change the gene id into gene name in go_inf

            return go_inf_filtered_geneName


        else:
            logger.info("get chartReport failed")
            raise Exception("upload genes to DAVID failed")


    def _checkSuccess(self,res):
        if res.find("DAVID: Functional Annotation Tools")==-1:
            return False
        else:
            return True

    '''
    pycurl helper class

    '''
    class PycurlHelper:

        def __init__(self):
            self.curl = pycurl.Curl()

        def get(self,url):
            buffer = StringIO()
            self.curl.setopt(pycurl.URL, url)
            self.curl.setopt(pycurl.COOKIEFILE, 'cookie.txt')
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
            self.curl.setopt(pycurl.CUSTOMREQUEST, "PUT")
            self.curl.setopt(self.curl.WRITEDATA, buffer)
            self.curl.setopt(pycurl.COOKIEJAR, 'cookie.txt')
            self.curl.setopt(pycurl.CAINFO, certifi.where())
            self.curl.perform()
            return buffer.getvalue().decode('iso-8859-1')

        def close(self):
            self.curl.close()



    def _uploadGene(self,pcHelper,idType,inputIds):

        data = [('idType', idType), ('uploadType', 'list'),('multiList','false'),('Mode','paste'),
                         ('useIndex','null'),('usePopIndex','null'),('demoIndex','null'),('ids',inputIds),
                         ('removeIndex','null'),('renameIndex','null'),('renamePopIndex','null'),('newName','null'),
                         ('combineIndex','null'),('selectedSpecies','null'),('uploadHTML','null'),('managerHTML','null'),
                         ('sublist',''),('rowids',''),('convertedListName','null'),('convertedPopName','null'),
                         ('pasteBox',inputIds),('Identifier',idType) , ('rbUploadType','list')]


        return pcHelper.sendMultipart(url="https://david.ncifcrf.gov/tools.jsp",data=data)




    def _changeGeneIdToNameInGO(self,go_inf_filtered,geneIdNameMapping):
        for i in range(0,len(go_inf_filtered)):

            geneNames = go_inf_filtered[i]["genes"].split(",")

            geneNameString = "|".join([geneIdNameMapping[geneName.strip().lower()] for geneName in geneNames]) 

            go_inf_filtered[i]["genes"] = geneNameString[:-1]# strip the last '|'

            return go_inf_filtered


    def _getGenesNamesByIds(self,pcHelper,geneIds):
        res = pcHelper.get('https://david.ncifcrf.gov/list.jsp')#get gene name and id
        parser = DavidDataScrawler.geneParser()
        parser.feed(res)#mapping between gene id and gene name

        genesIdName = parser.getParsedGeneNameWithId()
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

        return filterGO_inf





    # create a subclass and override the handler methods
    class GOParser(HTMLParser):
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



if __name__ == '__main__':
    app.run(debug="true",host="0.0.0.0")