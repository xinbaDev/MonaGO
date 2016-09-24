import unittest
from DavidDataScrawler import DavidDataScrawler
from DataProcess import DataProcess

class TestStringMethods(unittest.TestCase):
	def setUp(self):
		dataScrawler = DavidDataScrawler()
		DataProcess = DataProcess()

	def tearDown(self):
		return


	def test_Scrawler(self):
		davidScrawler.setParams(inputIds,idType,annotCat,pVal)
		
	


if __name__ == '__main__':
    unittest.main()