
import bgcArgoDMQC as bgc
import unittest

class configureTest(unittest.TestCase):

    def test_configure(self):

        bgc.configure.reset_config()
        bgc.configure.configure(
            operator='Christopher Gordon', 
            affiliation='Fisheries and Oceans Canada',
            orcid='0000-0002-1756-422X',
            argo_path='/Users/GordonC/Doduments/data/Argo/dac',
            ncep_path='/Users/GordonC/Doduments/data/NCEP',
            woa_path='/Users/GordonC/Doduments/data/WOA18',
        )
        config = bgc.configure.read_config()

        self.assertTrue(config['operator'] == 'Christopher Gordon')
        self.assertTrue(config['affiliation'] == 'Fisheries and Oceans Canada')

        bgc.configure.reset_config()
        config = bgc.configure.read_config()
        self.assertEqual(len(config.keys()), 0)

        bgc.configure.configure(
            operator='Christopher Gordon', 
            affiliation='Fisheries and Oceans Canada',
            orcid='0000-0002-1756-422X',
            argo_path='/Users/GordonC/Doduments/data/Argo/dac',
            ncep_path='/Users/GordonC/Doduments/data/NCEP',
            woa_path='/Users/GordonC/Doduments/data/WOA18',
        )

if __name__ == '__main__':
    unittest.main()
