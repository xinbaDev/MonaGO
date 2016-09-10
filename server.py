
#for pythonAPI
import sys
sys.path.append('../')
# import suds.metrics as metrics
# from suds import *
# from suds.client import Client
from datetime import datetime

#for praser
from HTMLParser import HTMLParser
import re

#for http request
import requests
# from requests_toolbelt.multipart.encoder import MultipartEncoder

from flask import Flask,render_template,request,send_from_directory
app = Flask(__name__)

#for pre processing the data
import dataProcessing

import json

import pycurl
from StringIO import StringIO
import certifi

import logging


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


with open("E:\\research\\zebrafish\\server\\js\\GO.js","r") as fr_GO:
	for i in fr_GO:
		GO_hier = json.loads(str(i)) 


with open('E:\\research\\zebrafish\\MonaGO\\MonaGO\\templates\\chord_layout.html',"r") as fr_html:
	html = "".join(fr_html.readlines())


metacount = 0
geneLists = {}
genesIdName = []
GO_hier_list = {}

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
		#html = davidPythonAPI(inputIds,idType,listName,listType,annotCat);

		annotCat = ','.join([annotCatDict[cat] for cat in annotCat.split(",")])

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

	with open("E:\\research\\zebrafish\\Data.txt","r") as fr:
		data = fr.readline()

	return data + html;


def getSessionId(s):
	r = s.get('http://david.abcc.ncifcrf.gov/tools.jsp')
	parser = sessionIdParser()
	parser.feed(r.text)
	return parser.returnSessionId()
	

def uploadGene(pcHelper,idType,inputIds):

	data = [('idType', idType), ('uploadType', 'list'),('multiList','false'),('Mode','paste'),
					 ('useIndex','null'),('usePopIndex','null'),('demoIndex','null'),('ids',inputIds),
					 ('removeIndex','null'),('renameIndex','null'),('renamePopIndex','null'),('newName','null'),
					 ('combineIndex','null'),('selectedSpecies','null'),('uploadHTML','null'),('managerHTML','null'),
					 ('sublist',''),('rowids',''),('convertedListName','null'),('convertedPopName','null'),
					 ('pasteBox',inputIds),('Identifier',idType) , ('rbUploadType','list')]


	return pcHelper.sendMultipart(url="https://david.ncifcrf.gov/tools.jsp",data=data)



def davidWebAPI(inputIds,idType,listName,listType,annotCat,pVal):
	global geneLists

	pcHelper = PycurlHelper()

	res = uploadGene(pcHelper,idType,inputIds)

	if _checkSuccess(res):
		url = 'https://david.ncifcrf.gov/chartReport.jsp?annot={0}&currentList=0'.format(annotCat)
		getGO_response = pcHelper.get(url)

		parser = GOParser()
		parser.feed(getGO_response)#get go_inf
		go_inf = parser.getGO_inf()



	else:
		logger.info("get chartReport failed")


	#sessionId = getSessionId(s)



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

	# go_hier = getGODependency(go_inf_reord)

	# #print go_hier

	# data = "<script>"+"var go_inf="+str(go_inf_reord)+";"+"var matrix="+matrix_count+";"+"var array_order="+array_order+";"+"var clusterHierData="+clusterHierData\
	# +";"+"var size="+str(len(go_inf_reord))+";"+"var goNodes="+str(go_hier)+"</script>"

	# # fw = open("E:\\research\\zebrafish\\Data.txt","w+")
	# # fw.write(data)
	# # fw.close()

	# return data+html

	return res;

def _checkSuccess(res):
	if res.find("DAVID: Functional Annotation Tools")==-1:
		return False
	else:
		return True


def filterGO(pVal,go_inf):

	filterGO_inf = [for i in go_inf if float(go_inf[i]["pVal"]) < float(pVal)]
	
	return filterGO_inf



# create a subclass and override the handler methods
class GOParser(HTMLParser):
	def __init__(self):
		self.go_inf = []

  # function to handle an opening tag in the doc
	def handle_starttag(self, tag, attrs):
		global geneLists
		global metacount

		#get GO_id,GO_name,p-value,count
		if tag == "a":

			m = re.search('(data/download/chart_\w+.txt)',attrs[0][1])
			if m!=None:
				url = 'http://david.abcc.ncifcrf.gov/'+m.group(0)
				s = requests.get(url)

				self._parseGO(s.text)

				

		#get gene rowid
		if tag == "img":
			if attrs[0][1] == 'graphics/two_tone_2_a.jpg':
				genes = attrs[6][1].split(";")[1]
				geneLists[metacount] = parseStringIntoList(genes)
				metacount+=1

	def _parseGO(self,GO):
		#get rid of first line
		line = GO.encode('ascii','ignore').split("\n")
		for i in range(1,len(line)-1):
			cell = line[i].split("\t")
			GO_term = cell[1].split("~")
			self.go_inf.append({"cat":cell[0],"GO_id":GO_term[0],"GO_name":GO_term[1],"count":cell[2],"pVal":cell[4],"genes":cell[5].lower().strip()})

	def getGO_inf(self):
		return self.go_inf

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






def parseStringIntoList(genes):
	index = genes.find("geneReport")
	geneList = genes[index+10:][2:-2].encode('ascii','ignore').split(",")
	return geneList



def getGODependency(GO_inf):
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