import sys
sys.path.append("..")

from unittest import TestCase, main
from scripts.parse_samplesheet import *
import yaml

class ParseSamplesheetTests(TestCase):
    samplesheets = None
    config = None

    def setUp(self):
        self.samplesheets = get_global_samplesheets('tests/data/')
        with open('../config.yaml', 'r') as f:
            self.config = yaml.load(f)

    def test_get_role(self):
        config = {'dirs': {'prefix': './',
                           'inputs': 'tests/',
                           'samplesheets': 'data'}}
        ukd_project = 'Alpss'
        with self.assertRaisesRegex(ValueError,
                                    ('Could not find an UKD project with name "%s". Available projects are:\n\t%s\n' % (ukd_project, '\n\t'.join(sorted(['AG_Remke', 'Alps', 'Fischer_Geron', 'Keimbahn', 'Maus_Hauer', 'ALL_family_LB']))))):
            get_role(ukd_project, None, None, self.samplesheets)

        ukd_project = 'Alps'
        ukd_entity_id = 'ALPS_6x'
        with self.assertRaisesRegex(ValueError,
                                    ('Could not find an UKD entity group with name "%s". Available entities for projects "%s" are:\n\t%s\n' % (ukd_entity_id, ukd_project, '\n\t'.join(sorted(['ALPS_64', 'ALPS_65', 'ALPS_66']))))):
            get_role(ukd_project, ukd_entity_id, None, self.samplesheets)

        ukd_entity_id = 'ALPS_66'
        ukd_entity_role = 'grandma'
        with self.assertRaisesRegex(ValueError,
                                    ('Could not find a role "%s" for UKD entity group with name "%s". Available roles are:\n\t%s\n' % (ukd_entity_role, ukd_entity_id, '\n\t'.join(sorted(['father', 'mother', 'patient']))))):
            get_role(ukd_project, ukd_entity_id, ukd_entity_role, self.samplesheets)

        ukd_entity_role = 'father'
        exp = 'Alps/ALPS_66_a'
        obs = get_role(ukd_project, ukd_entity_id, ukd_entity_role, self.samplesheets)
        self.assertEqual(exp, obs)

    def test_get_reference_exometrack(self):
        exp = {'Agilent_SureSelect_V5plusUTR.bed': [
            'Keimbahn/KB0198_c',
            'Keimbahn/KB0198_f',
            'Keimbahn/KB0198_m',
            'Keimbahn/KB0005_m',
            'Keimbahn/KB0043_c',
            'Keimbahn/KB0043_m',
            'Keimbahn/KB0077_c',
            'Keimbahn/KB0077_m',
            'Keimbahn/KB0122_c',
            'Keimbahn/KB0122_f',
            'Keimbahn/KB0122_m',
            'Alps/ALPS_64',
            'Alps/ALPS_64_a',
            'Alps/ALPS_64_b',
            'Alps/ALPS_65',
            'Alps/ALPS_65_a',
            'Alps/ALPS_65_b',
            'Keimbahn/KB0135_c',
            'Keimbahn/KB0135_f',
            'Keimbahn/KB0135_m',
            'Keimbahn/KB0136_c',
            'Keimbahn/KB0136_m',
            'Keimbahn/KB0199_c',
            'Keimbahn/KB0199_f',
            'Keimbahn/KB0199_m',
            'Keimbahn/KB0190_c',
            'Keimbahn/KB0190_m',
            'Keimbahn/KB0174_c',
            'Keimbahn/KB0174_m',
            'Keimbahn/KB0206_c',
            'Keimbahn/KB0206_f',
            'Keimbahn/KB0145_c',
            'Keimbahn/KB0145_m',
            'Alps/ALPS_66',
            'Alps/ALPS_66_a',
            'Alps/ALPS_66_b',
            'Fischer_Geron/F49',
            'Fischer_Geron/42'],
        'Agilent_mouseexome_S0276129_Regions_mm10_modified_merged.bed': [
            'Maus_Hauer/mAID_292_N',
            'Maus_Hauer/mAID_292_T',
            'Maus_Hauer/mAID_293_N',
            'Maus_Hauer/mAID_293_T',
            'Maus_Hauer/mMicro295_N',
            'Maus_Hauer/mMicro295_T',
            'Maus_Hauer/mMicro296_N',
            'Maus_Hauer/mMicro296_T',
            'Maus_Hauer/mMicro297_N',
            'Maus_Hauer/mMicro297_T',
            'Maus_Hauer/mMicro298_N',
            'Maus_Hauer/mMicro298_T',
            'Maus_Hauer/286C',
            'Maus_Hauer/286T',
            'Maus_Hauer/287C',
            'Maus_Hauer/287T',
            'Maus_Hauer/288C',
            'Maus_Hauer/288T',
            'Maus_Hauer/291C',
            'Maus_Hauer/291T'],
        "Agilent_SureSelect_Human_AllExon_V6_r2.hGRC37.Covered.bed": [
            'ALL_family_LB/LB_4_2',
            'ALL_family_LB/LB_4_3',
            'ALL_family_LB/LB_5_1',
            'ALL_family_LB/LB_5_2',
            'ALL_family_LB/LB_5_3']}

        for entity in self.samplesheets[~self.samplesheets['Sample_Project'].isin(['AG_Remke'])]['fastq-prefix'].unique():
            obs = get_reference_exometrack(entity, self.samplesheets, self.config)
            self.assertIn(entity, exp[obs])

    def test_get_bwa_mem_header(self):
        self.assertIn('SureSelectXTV5plusUTRautomated', get_bwa_mem_header('Alps/ALPS_64', self.samplesheets, self.config))
        self.assertIn('SureSelectXTmouse', get_bwa_mem_header('Maus_Hauer/mMicro297_T', self.samplesheets, self.config))
        self.assertIn('SureSelectHumanAllExonV6r2', get_bwa_mem_header('ALL_family_LB/LB_5_3', self.samplesheets, self.config))


if __name__ == '__main__':
    main()
