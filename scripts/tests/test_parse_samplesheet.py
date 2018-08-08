import sys
sys.path.append("..")

from unittest import TestCase, main
from scripts.parse_samplesheet import *

class ParseSamplesheetTests(TestCase):
    fp_samplesheet = 'tests/data/180614_SN737_0438_BCC7MCACXX_ukd.csv'

    def test_get_role(self):
        config = {'dirs': {'prefix': './',
                           'inputs': 'tests/',
                           'samplesheets': 'data'}}
        ukd_project = 'Alpss'
        with self.assertRaisesRegex(ValueError,
                                    ('Could not find an UKD project with name "%s". Available projects are:\n\t%s\n' % (ukd_project, '\n\t'.join(sorted(['AG_Remke', 'Alps', 'Fischer_Geron', 'Keimbahn', 'Maus_Hauer', 'Mouse_Spain']))))):
            get_role(ukd_project, None, None, config)

        ukd_project = 'Alps'
        ukd_entity_id = 'ALPS_6x'
        with self.assertRaisesRegex(ValueError,
                                    ('Could not find an UKD entity group with name "%s". Available entities for projects "%s" are:\n\t%s\n' % (ukd_entity_id, ukd_project, '\n\t'.join(sorted(['ALPS_64', 'ALPS_65', 'ALPS_66']))))):
            get_role(ukd_project, ukd_entity_id, None, config)

        ukd_entity_id = 'ALPS_66'
        ukd_entity_role = 'grandma'
        with self.assertRaisesRegex(ValueError,
                                    ('Could not find a role "%s" for UKD entity group with name "%s". Available roles are:\n\t%s\n' % (ukd_entity_role, ukd_entity_id, '\n\t'.join(sorted(['father', 'mother', 'patient']))))):
            get_role(ukd_project, ukd_entity_id, ukd_entity_role, config)

        ukd_entity_role = 'father'
        exp = '180614_SN737_0438_BCC7MCACXX/Alps/ALPS_66_a'
        obs = get_role(ukd_project, ukd_entity_id, ukd_entity_role, config)
        self.assertEqual(exp, obs)


if __name__ == '__main__':
    main()
