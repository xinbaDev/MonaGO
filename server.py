
import sys

from datetime import datetime

#for praser


from flask import Flask,render_template,request,send_from_directory


#for pre processing the data
from DataProcess import DataProcess
from DavidDataScrawler import DavidDataScrawler

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
    dataProcess = DataProcess()

    preProcessedData = dataProcess.dataProcess(go_inf_filtered_geneName)

    matrix_count = str(preProcessedData["matrix"]["matrix_count"])

    array_order = str(preProcessedData["gen_anno_reord"])

    clusterHierData = str(preProcessedData["clusterHierData"])

    go_inf_reord = preProcessedData["go_inf"]

    go_hier = preProcessedData["go_hier"]

    return matrix_count,array_order,go_hier,go_inf_reord,clusterHierData


if __name__ == '__main__':
    app.run(debug="true",host="0.0.0.0")