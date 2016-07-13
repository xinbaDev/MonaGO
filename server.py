
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

from flask import Flask,render_template,request,send_from_directory
app = Flask(__name__)

#for pre processing the data
import dataProcessing

import json

fr_GO = open("E:\\research\\zebrafish\\server\\js\\GO.js","r")

fr = open("E:\\research\\zebrafish\\Data.txt","r")

for i in fr_GO:
    GO_hier = json.loads(str(i))

fr_GO.close()

metacount = 0
geneLists = {}
go_inf=[];
genesIdName = []

@app.route('/')
def hello():
    return "hello"

@app.route('/test', methods=['POST','GET'])
def test():
    if request.method == 'GET':
        return render_template("test.html")
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

    return data + render_template("chord_layout.html")

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




def davidWebAPI(inputIds,idType,listName,listType,annotCat,pVal):
    global geneLists
    global go_inf

    s = requests.Session()

    get_rowsId_response = s.get("http://david.abcc.ncifcrf.gov/api.jsp?type="+idType+"&ids="+inputIds+"&tool=chartReport&annot="+annotCat);

    m = re.search("Request-URI Too Long",get_rowsId_response.text)

    if(m!=None):
        return "The requested URL's length exceeds the capacity"

    rowids = ""
    m = re.search('document.apiForm.rowids.value="(.+)"',get_rowsId_response.text)
    if(m!=None):
         rowids = m.group(1)
    m = re.search('document.apiForm.annot.value="(.+)"',get_rowsId_response.text)
    if(m!=None):
         annot = m.group(1)


    if(rowids==""):
        return "Less than 80% of your list has mapped to your chosen identifier type. Please use the Gene Conversion Tool on DAVID to determine the identifier type."

    # print 'http://david.abcc.ncifcrf.gov/chartReport.jsp?rowids='+rowids+"&annot="+annot

    getGO_response = s.get('http://david.abcc.ncifcrf.gov/chartReport.jsp?rowids='+rowids+"&annot="+annot)

    m = re.search("Request-URI Too Long",getGO_response.text)
    if(m!=None):
        return "The requested URL's length exceeds the capacity"

    parser = GOParser()

    parser.feed(getGO_response.text)#get go_inf

    go_inf = filterGO(pVal)

    rowids = set([])
    for i in range(0,len(geneLists)-1):
        for id in geneLists[i]:
            rowids.add(id)

    rowidstr = ""
    for i in rowids:
        rowidstr += (str(i)+",")#get gene rowid

    response = s.get('https://david.ncifcrf.gov/geneReport.jsp?rowids='+rowidstr)

    parser = geneParser()
    parser.feed(response.text)#mapping between gene id and gene name

    geneMapping = {};
    for i in range(0,len(genesIdName)-1,2):
        geneIds = genesIdName[i].split(",")
        if len(geneIds) > 1:
            for j in geneIds:
                geneMapping[j.strip().lower()] = genesIdName[i+1]
        else:
            geneMapping[genesIdName[i].lower()] = genesIdName[i+1]

    for i in range(0,len(go_inf)):
        genes = go_inf[i]["genes"].split(",")
        newString = ""
        for j in genes:
            newString += geneMapping[j.strip().lower()]+"|"
        go_inf[i]["genes"] = newString[:-1]#map gene id to gene name


    preProcessedData = dataProcessing.dataPreprocess(go_inf)

    matrix_count = str(preProcessedData["matrix"]["matrix_count"])

    array_order = str(preProcessedData["gen_anno_reord"])

    clusterHierData = str(preProcessedData["clusterHierData"])

    go_inf_reord = preProcessedData["go_inf"]

    go_hier = getGOdDependency(go_inf_reord)

    #print go_hier

    data = "<script>"+"var go_inf="+str(go_inf_reord)+";"+"var matrix="+matrix_count+";"+"var array_order="+array_order+";"+"var clusterHierData="+clusterHierData\
    +";"+"var size="+str(len(go_inf_reord))+";"+"var goNodes="+str(go_hier)+"</script>"

    fw = open("E:\\research\\zebrafish\\Data.txt","w+")
    fw.write(data)
    fw.close()

    return data+render_template("chord_layout.html")

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