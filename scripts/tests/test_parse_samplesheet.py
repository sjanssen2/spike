import sys
sys.path.append("..")

from unittest import TestCase, main
from scripts.parse_samplesheet import *

class ParseSamplesheetTests(TestCase):
    fp_samplesheet = 'tests/data/180614_SN737_0438_BCC7MCACXX_ukd.csv'

    def test_get_fastq_filenames(self):
        exp = ["Undetermined_L001_R1_001.fastq.gz",
               "Undetermined_L001_R2_001.fastq.gz",
               "AG_Remke/Chri_3/297_L001_R1_001.fastq.gz",
               "AG_Remke/Chri_3/297_L001_R2_001.fastq.gz",
               "AG_Remke/Chri_4/298_L001_R1_001.fastq.gz",
               "AG_Remke/Chri_4/298_L001_R2_001.fastq.gz",
               "AG_Remke/NPTH_17/299_L001_R1_001.fastq.gz",
               "AG_Remke/NPTH_17/299_L001_R2_001.fastq.gz",
               "AG_Remke/NPTH_18/300_L001_R1_001.fastq.gz",
               "AG_Remke/NPTH_18/300_L001_R2_001.fastq.gz",
               "AG_Remke/CLIS_10/277_L002_R1_001.fastq.gz",
               "AG_Remke/CLIS_10/277_L002_R2_001.fastq.gz",
               "AG_Remke/CLIS_17/278_L002_R1_001.fastq.gz",
               "AG_Remke/CLIS_17/278_L002_R2_001.fastq.gz",
               "AG_Remke/CLIS_8_b/279_L002_R1_001.fastq.gz",
               "AG_Remke/CLIS_8_b/279_L002_R2_001.fastq.gz",
               "AG_Remke/CLIS_13_b/280_L002_R1_001.fastq.gz",
               "AG_Remke/CLIS_13_b/280_L002_R2_001.fastq.gz",
               "AG_Remke/CLIS_21/272_L007_R1_001.fastq.gz",
               "AG_Remke/CLIS_21/272_L007_R2_001.fastq.gz",
               "AG_Remke/CLIS_22/273_L007_R1_001.fastq.gz",
               "AG_Remke/CLIS_22/273_L007_R2_001.fastq.gz",
               "AG_Remke/CLIS_23/274_L007_R1_001.fastq.gz",
               "AG_Remke/CLIS_23/274_L007_R2_001.fastq.gz",
               "AG_Remke/NPTH_1/281_L007_R1_001.fastq.gz",
               "AG_Remke/NPTH_1/281_L007_R2_001.fastq.gz",
               "AG_Remke/CLIS_1/263_L008_R1_001.fastq.gz",
               "AG_Remke/CLIS_1/263_L008_R2_001.fastq.gz",
               "AG_Remke/CLIS_2/264_L008_R1_001.fastq.gz",
               "AG_Remke/CLIS_2/264_L008_R2_001.fastq.gz",
               "AG_Remke/CLIS_5/265_L008_R1_001.fastq.gz",
               "AG_Remke/CLIS_5/265_L008_R2_001.fastq.gz",
               "AG_Remke/CLIS_24/275_L008_R1_001.fastq.gz",
               "AG_Remke/CLIS_24/275_L008_R2_001.fastq.gz",
               "Undetermined_L002_R1_001.fastq.gz",
               "Undetermined_L002_R2_001.fastq.gz",
               "Undetermined_L003_R1_001.fastq.gz",
               "Undetermined_L003_R2_001.fastq.gz",
               "Undetermined_L004_R1_001.fastq.gz",
               "Undetermined_L004_R2_001.fastq.gz",
               "Undetermined_L005_R1_001.fastq.gz",
               "Undetermined_L005_R2_001.fastq.gz",
               "Alps/ALPS_66_L005_R1_001.fastq.gz",
               "Alps/ALPS_66_L005_R2_001.fastq.gz",
               "Alps/ALPS_66_a_L005_R1_001.fastq.gz",
               "Alps/ALPS_66_a_L005_R2_001.fastq.gz",
               "Alps/ALPS_66_b_L005_R1_001.fastq.gz",
               "Alps/ALPS_66_b_L005_R2_001.fastq.gz",
               "Alps/ALPS_66_L006_R1_001.fastq.gz",
               "Alps/ALPS_66_L006_R2_001.fastq.gz",
               "Alps/ALPS_66_a_L006_R1_001.fastq.gz",
               "Alps/ALPS_66_a_L006_R2_001.fastq.gz",
               "Alps/ALPS_66_b_L006_R1_001.fastq.gz",
               "Alps/ALPS_66_b_L006_R2_001.fastq.gz",
               "Fischer_Geron/F49_L005_R1_001.fastq.gz",
               "Fischer_Geron/F49_L005_R2_001.fastq.gz",
               "Fischer_Geron/42_L005_R1_001.fastq.gz",
               "Fischer_Geron/42_L005_R2_001.fastq.gz",
               "Fischer_Geron/F49_L006_R1_001.fastq.gz",
               "Fischer_Geron/F49_L006_R2_001.fastq.gz",
               "Fischer_Geron/42_L006_R1_001.fastq.gz",
               "Fischer_Geron/42_L006_R2_001.fastq.gz",
               "Undetermined_L006_R1_001.fastq.gz",
               "Undetermined_L006_R2_001.fastq.gz",
               "Undetermined_L007_R1_001.fastq.gz",
               "Undetermined_L007_R2_001.fastq.gz",
               "Undetermined_L008_R1_001.fastq.gz",
               "Undetermined_L008_R2_001.fastq.gz",
               'Maus_Hauer/286C_L003_R1_001.fastq.gz',
               'Maus_Hauer/286C_L003_R2_001.fastq.gz',
               'Maus_Hauer/286T_L003_R1_001.fastq.gz',
               'Maus_Hauer/286T_L003_R2_001.fastq.gz',
               'Maus_Hauer/287C_L003_R1_001.fastq.gz',
               'Maus_Hauer/287C_L003_R2_001.fastq.gz',
               'Maus_Hauer/287T_L003_R1_001.fastq.gz',
               'Maus_Hauer/287T_L003_R2_001.fastq.gz',
               'Maus_Hauer/288C_L004_R1_001.fastq.gz',
               'Maus_Hauer/288C_L004_R2_001.fastq.gz',
               'Maus_Hauer/288T_L004_R1_001.fastq.gz',
               'Maus_Hauer/288T_L004_R2_001.fastq.gz',
               'Maus_Hauer/291C_L004_R1_001.fastq.gz',
               'Maus_Hauer/291C_L004_R2_001.fastq.gz',
               'Maus_Hauer/291T_L004_R1_001.fastq.gz',
               'Maus_Hauer/291T_L004_R2_001.fastq.gz']

        obs = get_fastq_filenames(self.fp_samplesheet)
        self.assertEqual(sorted(exp), sorted(obs))

    def test_get_sample_fastqprefixes(self):
        exp = [
            'AG_Remke/CLIS_1/263',
            'Alps/ALPS_66',
            'Maus_Hauer/291T',
            'Maus_Hauer/286T',
            'AG_Remke/CLIS_13_b/280',
            'Alps/ALPS_66_a',
            'AG_Remke/CLIS_21/272',
            'Maus_Hauer/288T',
            'Fischer_Geron/F49',
            'AG_Remke/CLIS_5/265',
            'Maus_Hauer/287C',
            'AG_Remke/Chri_3/297',
            'AG_Remke/CLIS_22/273',
            'AG_Remke/CLIS_8_b/279',
            'AG_Remke/CLIS_17/278',
            'Maus_Hauer/288C',
            'Fischer_Geron/42',
            'AG_Remke/CLIS_10/277',
            'AG_Remke/Chri_4/298',
            'AG_Remke/NPTH_1/281',
            'AG_Remke/CLIS_23/274',
            'AG_Remke/CLIS_2/264',
            'Alps/ALPS_66_b',
            'AG_Remke/NPTH_18/300',
            'Maus_Hauer/287T',
            'AG_Remke/CLIS_24/275',
            'AG_Remke/NPTH_17/299',
            'Maus_Hauer/291C',
            'Maus_Hauer/286C']
        obs = get_sample_fastqprefixes(self.fp_samplesheet)
        self.assertEqual(sorted(exp), sorted(obs))

    # def test_get_lanes_for_sampleID(self):
    #     data = [
    #         ('Chri_3', '297', '1', [1]),
    #         ('Chri_4', '298', '2', [1]),
    #         ('NPTH_17', '299', '3', [1]),
    #         ('NPTH_18', '300', '4', [1]),
    #         ('CLIS_10', '277', '5', [2]),
    #         ('CLIS_17', '278', '6', [2]),
    #         ('CLIS_8_b', '279', '7', [2]),
    #         ('CLIS_13_b', '280', '8', [2]),
    #         ('286C', '', '9', [3]),
    #         ('286T', '', '10', [3]),
    #         ('287C', '', '11', [3]),
    #         ('287T', '', '12', [3]),
    #         ('288C', '', '13', [4]),
    #         ('288T', '', '14', [4]),
    #         ('291C', '', '15', [4]),
    #         ('291T', '', '16', [4]),
    #         ('ALPS_66', '', '17', [5, 6]),
    #         ('ALPS_66_a', '', '18', [5, 6]),
    #         ('ALPS_66_b', '', '19', [5, 6]),
    #         ('F49', '', '20', [5, 6]),
    #         ('42', '', '21', [5, 6]),
    #         ('CLIS_21', '272', '22', [7]),
    #         ('CLIS_22', '273', '23', [7]),
    #         ('CLIS_23', '274', '24', [7]),
    #         ('NPTH_1', '281', '25', [7]),
    #         ('CLIS_1', '263', '26', [8]),
    #         ('CLIS_2', '264', '27', [8]),
    #         ('CLIS_5', '265', '28', [8]),
    #         ('CLIS_24', '275', '29', [8])]
    #
    #     for (sampleName, sampleID, sidx, exp) in data:
    #         obs = get_lanes_for_sampleID(self.fp_samplesheet, sampleName + ('/' if sampleID != '' else ''), sampleID, sidx)
    #         self.assertEqual(sorted(exp), sorted(obs))

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
