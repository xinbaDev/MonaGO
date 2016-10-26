import pycurl
from StringIO import StringIO
import certifi
import logging
import time
import select

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