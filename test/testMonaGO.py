import unittest, sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'server'))
from DavidDataScrawler import DavidDataScrawler
# from DataProcess import DataProcess


class Test(unittest.TestCase):
	def setUp(self):
		self.davidScrawler = DavidDataScrawler()
		# DataProcess = DataProcess()
		

	def tearDown(self):
		pass

	def test_Scrawler_success(self):
		inputIds = "1007_s_at,1053_at,117_at,121_at,1255_g_at,1294_at,1316_at,1320_at,1405_i_at,1431_at,1438_at,1487_at,1494_f_at,1598_g_at"
		idType = "AFFYMETRIX_3PRIME_IVT_ID"
		annotCat = "25"
		pVal = "0.05"
		self.davidScrawler.setParams(inputIds,idType,annotCat,pVal)
		go = self.davidScrawler.run()
		self.assertTrue(len(go) > 0)

	def test_Scrawler_inputIds_isEmpty(self):
		inputIds=""
		idType = "AFFYMETRIX_3PRIME_IVT_ID"
		annotCat = "25"
		pVal = "0.05"
		self.davidScrawler.setParams(inputIds,idType,annotCat,pVal)
		with self.assertRaises(Exception):
			go = self.davidScrawler.run()

	


if __name__ == '__main__':
    unittest.main()