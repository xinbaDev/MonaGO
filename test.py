import unittest
from DavidDataScrawler import DavidDataScrawler
from twisted.internet import defer

class Test(unittest.TestCase):
    def setUp(self):
        self.davidScrawler = DavidDataScrawler()
        # DataProcess = DataProcess()
        self.inputIds = "38926_at,38236_at,35367_at,1910_s_at,1391_s_at,38674_at,31576_at,606_at,34703_f_at"
        self.idType = "AFFYMETRIX_3PRIME_IVT_ID"
        self.annotCat = "25"
        self.pVal = "0.05"

    def tearDown(self):
        return

    @defer.inlineCallbacks
    def test_Scrawler_Success(self):
        self.davidScrawler.setParams(self.inputIds,self.idType,self.annotCat,self.pVal)
        go = yield self.davidScrawler.run()
        self.assertTrue(go)#go is not empty

    @defer.inlineCallbacks
    def test_Scrawler_inputIds_isEmpty(self):
        self.inputIds=""
        self.davidScrawler.setParams(self.inputIds,self.idType,self.annotCat,self.pVal)
        with self.assertRaises(Exception):
            yield self.davidScrawler.run()

    @defer.inlineCallbacks
    def test_Scrawler_wrong_idType(self):
        self.idType = "General Symbol"
        self.davidScrawler.setParams(self.inputIds,self.idType,self.annotCat,self.pVal)
        with self.assertRaises(Exception):
            yield self.davidScrawler.run()

    @defer.inlineCallbacks
    def test_Scrawler_wrong_annotCat(self):
        self.annotCat = "99999" #inexistent annoation type
        self.davidScrawler.setParams(self.inputIds,self.idType,self.annotCat,self.pVal)
        with self.assertRaises(Exception):
            yield self.davidScrawler.run()

    @defer.inlineCallbacks
    def test_Scrawler_small_pVal(self):
        self.pVal = "0.00005"
        self.davidScrawler.setParams(self.inputIds,self.idType,self.annotCat,self.pVal)
        with self.assertRaises(Exception):
            yield self.davidScrawler.run()


if __name__ == '__main__':
    unittest.main()