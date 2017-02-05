import sys

from datetime import datetime

#for praser
from flask import Flask,render_template,request,send_from_directory,Response

#for pre processing the data
from DataProcess import DataProcess
from DataProcess2 import DataProcess2
from DavidDataScrawler import DavidDataScrawler

import logging
import time

import base64

from logTime import logTime

import ast

logging.basicConfig(filename="debug.txt",level=logging.INFO)
logger = logging.getLogger(__name__)
remote_server = False;

if(remote_server):
    root_dir = "/root/alex/myproject/"
else:
    root_dir = ""

app = Flask(__name__)
@app.route('/', methods=['POST','GET'])
def index():

    annotCatDict = {
        'GOTERM_BP_FAT':'25',#biological process
        'GOTERM_CC_FAT':'32',#celluar component
        'GOTERM_MF_FAT':'39'#melocular function
    }

    if request.method == 'GET':
        return render_template("index.html")
    else:
        if request.form['inputGOs']:
            go = ast.literal_eval(request.form['inputGOs'])

        else:
            #parameters needed for querying DAVID
            inputIds = request.form['inputIds']
            idType = request.form['idType']
            annotCat = request.form['annotCat']
            pVal = request.form['pVal'];
            #transform annotation name to number recognized by DAVID(e.g. GOTERM_BP_FAT to 25) .
            annotCat = ','.join([annotCatDict[cat] for cat in annotCat.split(",")])

            go,status = getDataFromDavid(inputIds,idType,annotCat,pVal)

            if status == False:
                return "Failure to get data, please make sure the identifier is correct"

        matrix_count,array_order,go_hier,go_inf_reord,clusterHierData = processedData(go)

        if not matrix_count:
            return "Failure to process data"

        with open(root_dir+'templates/chord_layout.html',"r") as fr_html:
            html = "".join(fr_html.readlines())

        data = "<script>"+"var go_inf="+str(go_inf_reord)+";"+"var matrix="+str(matrix_count)+";"+"var array_order="+str(array_order)+";"\
        +"var clusterHierData="+str(clusterHierData) +";"+"var size="+str(len(go_inf_reord))+";"+"var goNodes="+str(go_hier)+"</script>"

        return data+html

@app.route('/index2', methods=['POST','GET'])
def index2():

    annotCatDict = {
        'GOTERM_BP_FAT':'25',#biological process
        'GOTERM_CC_FAT':'32',#celluar component
        'GOTERM_MF_FAT':'39'#melocular function
    }

    if request.method == 'GET':
        return render_template("index2.html")
    else:

        if request.form['inputGOs']:
            go = request.form['inputGOs']

        else:
        #parameters needed for querying DAVID
            inputIds = request.form['inputIds']
            idType = request.form['idType']
            annotCat = request.form['annotCat']
            pVal = request.form['pVal'];
            #transform annotation name to number recognized by DAVID(e.g. GOTERM_BP_FAT to 25) .
            annotCat = ','.join([annotCatDict[cat] for cat in annotCat.split(",")])

            go,status = getDataFromDavid(inputIds,idType,annotCat,pVal)

            if status == False:
                return "Failure to get data, please make sure the identifier is correct"

        matrix_count,array_order,go_hier,go_inf_reord,clusterHierData = processedData2(go)

        if not matrix_count:
            return "Failure to process data"

        with open(root_dir+'templates/chord_layout.html',"r") as fr_html:
            html = "".join(fr_html.readlines())

        data = "<script>"+"var go_inf="+str(go_inf_reord)+";"+"var matrix="+str(matrix_count)+";"+"var array_order="+str(array_order)+";"\
        +"var clusterHierData="+str(clusterHierData) +";"+"var size="+str(len(go_inf_reord))+";"+"var goNodes="+str(go_hier)+"</script>"

        return data+html

@app.route('/css/<fileName>')
def getCss(fileName):
    return send_from_directory(root_dir+'css', fileName)

@app.route('/css/font-awesome-4.6.3/css/<fileName>')
def getFontAwesome(fileName):
    return send_from_directory(root_dir+'css/font-awesome-4.6.3/css/', fileName)

@app.route('/css/font-awesome-4.6.3/fonts/<fileName>')
def getFont(fileName):
    return send_from_directory(root_dir+'css/font-awesome-4.6.3/fonts/', fileName)

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

    return data + html

@app.route('/my/img/my_logo.jpg')
def getMyLogo():
    with open("visitors.txt",'a') as fw:
        fw.write("remote address: {}  time: {}\n".format(request.remote_addr,datetime.today()))
    return send_from_directory(root_dir+'my/img','my_logo.jpg')

@app.route('/help.html')
def getHelp():
    return send_from_directory(root_dir+'templates','help.html')

@app.route('/getPic',methods=['POST','GET'])
def getPic():
    b64_string = request.form['png'].split(",")[1]

    b64_string += "=" * ((4 - len(request.form['png']) % 4) % 4)
    return Response(
        base64.b64decode(b64_string),
        mimetype="image/png",
        headers={"Content-disposition":
                 "attachment; filename=chart.png"})

@app.route('/export',methods=['POST','GET'])
def export():
    string = request.form['file']

    return Response(
        string,
        mimetype="file/txt",
        headers={"Content-disposition":
                 "attachment; filename=export.monago"})



@logTime
def getDataFromDavid(inputIds,idType,annotCat,pVal):
    '''
    send https request to David and get GO information
    
    Args:
        inputIds:a list of gene ids
        idType:the type of gene ids, such as AFFYMETRIX_3PRIME_IVT_ID
        annotcat: type of annotation wanted, such as GO
        pVal: p-value

    Return:
        A list of GO terms
    '''
    davidScrawler = DavidDataScrawler()
    davidScrawler.setParams(inputIds,idType,annotCat,pVal)

    try:
        go = davidScrawler.run()
    except Exception as e:
        logger.error(str(e))
        return [],False
    else:
        return go,True


@logTime
def processedData(go):
    '''
    generate necessary data for visualization

    Args:
        a list of GO terms

    Return:
        matrix:a matrix M where M(i,j) represents the number of intersected genes betweeen GO[i] and GO[j]
        go_index_reord:an array representing the position change of GO terms after hieracical clustering
        go_hier:a list of GO that are ancesters of enriched GO terms.
        go_inf_reord:an array of enriched GO terms
        clusterHierData:an array storing hierarcical data use to generated hierarcical tree
    '''
    with open("datafromdavid.txt","w") as fw:
        fw.write(str(go))

    dataProcess = DataProcess()
    try:
        preProcessedData = dataProcess.dataProcess(go)

    except Exception as e:
        logger.error(str(e))
    else:
        matrix = preProcessedData["matrix"]["matrix_count"]
        go_index_reord = preProcessedData["go_index_reord"]
        go_hier = preProcessedData["go_hier"]
        go_inf_reord = preProcessedData["go_inf"]
        clusterHierData = preProcessedData["clusterHierData"]
        
        return matrix,go_index_reord,go_hier,go_inf_reord,clusterHierData

@logTime
def processedData2(go):
    '''
    generate necessary data for visualization

    Args:
        a list of GO terms

    Return:
        matrix:a matrix M where M(i,j) represents the number of intersected genes betweeen GO[i] and GO[j]
        go_index_reord:an array representing the position change of GO terms after hieracical clustering
        go_hier:a list of GO that are ancesters of enriched GO terms.
        go_inf_reord:an array of enriched GO terms
        clusterHierData:an array storing hierarcical data use to generated hierarcical tree
    '''
    dataProcess = DataProcess2()

    try:
        preProcessedData = dataProcess.dataProcess(go)

    except Exception as e:
        logger.error(str(e))
    else:
        matrix = preProcessedData["matrix"]["matrix_count"]
        go_index_reord = preProcessedData["go_index_reord"]
        go_hier = preProcessedData["go_hier"]
        go_inf_reord = preProcessedData["go_inf"]
        clusterHierData = preProcessedData["clusterHierData"]
        
        return matrix,go_index_reord,go_hier,go_inf_reord,clusterHierData

if __name__ == '__main__':
    app.run(debug="true",host="0.0.0.0",threaded = True)