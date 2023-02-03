
import bgcArgoDMQC as bgc
import unittest

class configureTest(unittest.TestCase):

    def test_configure(self):

        bgc.configure.reset_config()
        bgc.configure.configure(operator='Christopher Gordon', affiliation='DFO')
        config = bgc.configure.read_config()

        self.assertTrue(config['operator'] == 'Christopher Gordon')
        self.assertTrue(config['affiliation'] == 'DFO')

        bgc.configure.reset_config()
        config = bgc.configure.read_config()
        self.assertEqual(len(config.keys() == 0))

if __name__ == '__main__':
    unittest.main()
