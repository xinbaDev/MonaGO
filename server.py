
#for pythonAPI
import sys
sys.path.append('../')
import suds.metrics as metrics
from suds import *
from suds.client import Client
from datetime import datetime

#for praser
from HTMLParser import HTMLParser
import re

#for http request
import requests
from requests_toolbelt.multipart.encoder import MultipartEncoder

from flask import Flask,render_template,request,send_from_directory
app = Flask(__name__)

#for pre processing the data
import dataProcessing

import json

import pycurl
from StringIO import StringIO
import certifi

fr_GO = open("E:\\research\\zebrafish\\server\\js\\GO.js","r")

fr = open("E:\\research\\zebrafish\\Data.txt","r")
fr_html = open('E:\\research\\zebrafish\\MonaGO\\MonaGO\\templates\\chord_layout.html',"r")

for i in fr_GO:
    GO_hier = json.loads(str(i))

html = "".join(fr_html.readlines())

fr_GO.close()
fr_html.close()

metacount = 0
geneLists = {}
go_inf=[];
genesIdName = []

@app.route('/')
def hello():
    return "hello"

@app.route('/index', methods=['POST','GET'])
def index():
    if request.method == 'GET':
        return render_template("index.html")
    else:

        inputIds = request.form['inputIds']
        idType = request.form['idType']
        listName = "test"
        listType = request.form['listType']
        annotCat = request.form['annotCat']
        pVal = request.form['pVal'];
        #html = davidPythonAPI(inputIds,idType,listName,listType,annotCat);

        html = davidWebAPI(inputIds,idType,listName,listType,annotCat,pVal);

        return html

@app.route('/test', methods='POST')
def test():
    request.get_data()
    print(request.data)

@app.route('/css/<fileName>')
def getCss(fileName):
    return send_from_directory('css', fileName)

@app.route('/img/<fileName>')
def getImg(fileName):
    return send_from_directory('img',fileName)

@app.route('/js/<fileName>')
def getJs(fileName):
    return send_from_directory('js', fileName)

@app.route('/fonts/<fileName>')
def getFonts(fileName):
    return send_from_directory('fonts', fileName)

@app.route('/txt/<fileName>')
def getText(fileName):
    return send_from_directory('text', fileName)

@app.route('/demo')
def returnDemo():
    
    data = fr.readline()

    return data + html;

def davidPythonAPI(inputIds,idType,listName,listType,annotCat):
        url = 'http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl'
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
        client.service.addList(inputIds, idType, listName, listType)

        #print client.service.getDefaultCategoryNames()

        client.service.setCategories(annotCat)

        return parsePythonAPIResponse(client.service.getChartReport(0.1,2))

def parsePythonAPIResponse(response):
    return str(response)


def getSessionId(s):
    r = s.get('http://david.abcc.ncifcrf.gov/tools.jsp')
    parser = sessionIdParser()
    parser.feed(r.text)
    return parser.returnSessionId()
    

def uploadGene(s,idType,inputIds):
    print inputIds
    print idType

    # payload = MultipartEncoder(
    #             [
    #                 ('idType', idType), ('uploadType', 'list'),('multiList','false'),('Mode','paste'),
    #                 ('useIndex','null'),('usePopIndex','null'),('demoIndex','null'),('ids',inputIds),
    #                 ('removeIndex','null'),('renameIndex','null'),('renamePopIndex','null'),('newName','null'),
    #                 ('combineIndex','null'),('selectedSpecies','null'),('uploadHTML','null'),('managerHTML','null'),
    #                 ('sublist',''),('rowids',''),('convertedListName','null'),('convertedPopName','null'),
    #                 ('pasteBox',inputIds),('Identifier',idType) , ('rbUploadType','list')
    #             ], "----WebKitFormBoundaryMghM0fh74hiYqk68"
                
    # )

    # print payload._calculate_length()

    # header = {
    #     'Accept':'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
    #     'Accept-Encoding':'gzip, deflate',
    #     'Accept-Language':'zh-CN,zh;q=0.8,en;q=0.6,zh-TW;q=0.4,ja;q=0.2,en-AU;q=0.2',
    #     'Cache-Control':'no-cache',
    #     'Connection':'keep-alive',
    #     'Pragma':'no-cache',
    #     'User-Agent':'Mozilla/5.0 (Windows NT 6.1; WOW64)',
    #     'Content-Type': payload.content_type,
    #     'Referer':'https://david.ncifcrf.gov/tools.jsp',
    #     'Host':'david.ncifcrf.gov',
    #     'Upgrade-Insecure-Requests':'1',
    #     'Origin':'https://david.ncifcrf.gov'
    # }

    myCookies = {
   
    }

    # data=dict(idType=idType,uploadType="list",multiList="false",Mode="paste",useIndex="null",usePopIndex="null",demoIndex='null',
    #                 ids=inputIds,removeIndex='null',renameIndex='null',renamePopIndex='null',newName='null',
    #                 combineIndex='null',selectedSpecies='null',uploadHTML='null',managerHTML='null',
    #                 sublist='',rowids='',convertedListName='null',convertedPopName='null',
    #                 pasteBox=inputIds,Identifier=idType, rbUploadType='list')

    data = [('idType', idType), ('uploadType', 'list'),('multiList','false'),('Mode','paste'),
                     ('useIndex','null'),('usePopIndex','null'),('demoIndex','null'),('ids',inputIds),
                     ('removeIndex','null'),('renameIndex','null'),('renamePopIndex','null'),('newName','null'),
                     ('combineIndex','null'),('selectedSpecies','null'),('uploadHTML','null'),('managerHTML','null'),
                     ('sublist',''),('rowids',''),('convertedListName','null'),('convertedPopName','null'),
                     ('pasteBox',inputIds),('Identifier',idType) , ('rbUploadType','list')]



    buffer = StringIO()
    c = pycurl.Curl()
    c.setopt(pycurl.URL, "https://david.ncifcrf.gov/tools.jsp")
    c.setopt(pycurl.HTTPPOST, data)
    c.setopt(pycurl.CUSTOMREQUEST, "PUT")
    c.setopt(c.WRITEDATA, buffer)
    c.setopt(pycurl.CAINFO, certifi.where())
    c.perform()
    c.close()
    #myCookies.update(s.cookies.get_dict())

    # r = s.post("http://david.abcc.ncifcrf.gov/tools.jsp", files=files)
    # myCookies.update(s.cookies.get_dict())
    
    body = buffer.getvalue()

    with open("E:\\research\\zebrafish\\result.html","w+") as f:
        f.write(body.decode('iso-8859-1'))


    return "test"


def davidWebAPI(inputIds,idType,listName,listType,annotCat,pVal):
    global geneLists
    global go_inf

    s = requests.Session()

    #get is not feasible
    #get_rowsId_response = s.get("http://david.abcc.ncifcrf.gov/api.jsp?type="+idType+"&ids="+inputIds+"&tool=chartReport&annot="+annotCat)
    # m = re.search("Request-URI Too Long",get_rowsId_response.text)

    # if(m!=None):
    #     return "The requested URL's length exceeds the capacity"

    # rowids = ""
    # m = re.search('document.apiForm.rowids.value="(.+)"',get_rowsId_response.text)
    # if(m!=None):
    #      rowids = m.group(1)
    # m = re.search('document.apiForm.annot.value="(.+)"',get_rowsId_response.text)
    # if(m!=None):
    #      annot = m.group(1)


    # if(rowids==""):
    #     return "Less than 80% of your list has mapped to your chosen identifier type. Please use the Gene Conversion Tool on DAVID to determine the identifier type."
        # print 'http://david.abcc.ncifcrf.gov/chartReport.jsp?rowids='+rowids+"&annot="+annot

    # getGO_response = s.get('http://david.abcc.ncifcrf.gov/chartReport.jsp?rowids='+rowids+"&annot="+annot)

    # m = re.search("Request-URI Too Long",getGO_response.text)
    # if(m!=None):
    #     return "The requested URL's length exceeds the capacity"

    #sessionId = getSessionId(s)

    res = uploadGene(s,idType,inputIds)

    fw = open("E:\\research\\zebrafish\\Data.html","wb")
    fw.write(res.encode('utf-8','strict'))
    fw.close()

    # getGO_response = s.get('https://david.ncifcrf.gov/chartReport.jsp?annot='+annotCat+'&currentList=0')
    # print getGO_response.text

    # parser = GOParser()

    # parser.feed(getGO_response.text)#get go_inf

    # #go_inf = filterGO(pVal)

    # # rowids = set([])
    # # for i in range(0,len(geneLists)-1):
    # #     for id in geneLists[i]:
    # #         rowids.add(id)

    # # rowidstr = ""
    # # for i in rowids:
    # #     rowidstr += (str(i)+",")#get gene rowid

    # # response = s.get('https://david.ncifcrf.gov/geneReport.jsp?rowids='+rowidstr)

    # geneName = s.get('https://david.ncifcrf.gov/list.jsp')#get gene name and id

    # parser = geneParser()
    # parser.feed(geneName.text)#mapping between gene id and gene name

    # geneMapping = {}
    # for i in range(0,len(genesIdName)-1,2):
    #     geneIds = genesIdName[i].split(",")
    #     if len(geneIds) > 1:
    #         for j in geneIds:
    #             geneMapping[j.strip().lower()] = genesIdName[i+1]
    #     else:
    #         geneMapping[genesIdName[i].lower()] = genesIdName[i+1]

    # for i in range(0,len(go_inf)):
    #     genes = go_inf[i]["genes"].split(",")
    #     newString = ""
    #     for j in genes:
    #         newString += geneMapping[j.strip().lower()]+"|"
    #     go_inf[i]["genes"] = newString[:-1]#map gene id to gene name

    # print go_inf

    # preProcessedData = dataProcessing.dataPreprocess(go_inf)

    # matrix_count = str(preProcessedData["matrix"]["matrix_count"])

    # array_order = str(preProcessedData["gen_anno_reord"])

    # clusterHierData = str(preProcessedData["clusterHierData"])

    # go_inf_reord = preProcessedData["go_inf"]

    # go_hier = getGOdDependency(go_inf_reord)

    # #print go_hier

    # data = "<script>"+"var go_inf="+str(go_inf_reord)+";"+"var matrix="+matrix_count+";"+"var array_order="+array_order+";"+"var clusterHierData="+clusterHierData\
    # +";"+"var size="+str(len(go_inf_reord))+";"+"var goNodes="+str(go_hier)+"</script>"

    # # fw = open("E:\\research\\zebrafish\\Data.txt","w+")
    # # fw.write(data)
    # # fw.close()

    # return data+html
    return "Mock";

def filterGO(pVal):
    filterGO_inf = []
    for i in range(0,len(go_inf)):

        if float(go_inf[i]["pVal"]) < float(pVal):
            filterGO_inf.append(go_inf[i])

    return filterGO_inf



# create a subclass and override the handler methods
class GOParser(HTMLParser):
  # function to handle an opening tag in the doc
    def handle_starttag(self, tag, attrs):
        global geneLists
        global metacount
        global go_inf

        #get GO_id,GO_name,p-value,count
        if tag == "a":

            m = re.search('(data/download/chart_\w+.txt)',attrs[0][1])
            if m!=None:
                url = 'http://david.abcc.ncifcrf.gov/'+m.group(0)
                s = requests.get(url)

                go_inf = parseGO(s.text)

                

        #get gene rowid
        if tag == "img":
            if attrs[0][1] == 'graphics/two_tone_2_a.jpg':
                genes = attrs[6][1].split(";")[1]
                geneLists[metacount] = parseStringIntoList(genes)
                metacount+=1




class geneParser(HTMLParser):
    table = 0
    tr = 0
    td = 0
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
        global genesIdName
        if self.tr == 1:
            if data != "RG" and self.td!= 4 and "\n" not in data:
                genesIdName.append(data.encode('ascii','ignore'))

class sessionIdParser(HTMLParser):
    def __init__(self):
        HTMLParser.__init__(self)
        self.sessionId = 0

    def handle_starttag(self, tag, attrs):
        if tag == 'input':
            if attrs[1][0] == 'name' and attrs[1][1] == 'SESSIONID':
                self.sessionId = attrs[2][1]

    def returnSessionId(self):
        return self.sessionId



def parseGO(GO):
    GO_array=[]
    #get rid of first line
    line = GO.encode('ascii','ignore').split("\n")
    for i in range(1,len(line)-1):
        cell = line[i].split("\t")
        GO_term = cell[1].split("~")
        GO_array.append({"cat":cell[0],"GO_id":GO_term[0],"GO_name":GO_term[1],"count":cell[2],"pVal":cell[4],"genes":cell[5].lower().strip()})
    return GO_array


def parseStringIntoList(genes):
    index = genes.find("geneReport")
    geneList = genes[index+10:][2:-2].encode('ascii','ignore').split(",")
    return geneList

GO_hier_list = {}

def getGOdDependency(GO_inf):
    global GO_hier_list
    GO_hier_list = {}
    for gos in GO_inf:
        recuriveGetGOId(gos["GO_id"].encode('ascii','ignore'))
    return json.dumps(GO_hier_list)


def recuriveGetGOId(GO_id):
    global GO_hier_list

    GO_hier_list[GO_id]=GO_hier[GO_id]

    for i in GO_hier[GO_id]["p"]:
        if not GO_hier_list.has_key(i):
            recuriveGetGOId(i.encode('ascii','ignore'))

if __name__ == '__main__':
    app.run(debug = True)