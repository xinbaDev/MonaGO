import unittest
from DavidDataScrawler import DavidDataScrawler
from DataProcess import DataProcess

class Test(unittest.TestCase):
	def setUp(self):
		dataScrawler = DavidDataScrawler()
		DataProcess = DataProcess()
		inputIds = "38926_at,38236_at,35367_at,1910_s_at,1391_s_at,38674_at,31576_at,606_at,34703_f_at"
		idType = "AFFYMETRIX_3PRIME_IVT_ID"
		annotCat = "25"
		pVal = "0.05"

	def tearDown(self):
		return



	def test_Scrawler_inputIds_isEmpty(self):
		inputIds=""
		davidScrawler.setParams(inputIds,idType,annotCat,pVal)
		
	


if __name__ == '__main__':
    unittest.main()