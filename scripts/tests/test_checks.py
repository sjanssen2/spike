import sys
sys.path.append("..")

from unittest import TestCase, main
import tempfile
import os
from scripts.checks import *

class ChecksTests(TestCase):
    fp_tmpfile = None

    def setUp(self):
        _, self.fp_tmpfile = tempfile.mkstemp()

    def tearDown(self):
        os.remove(self.fp_tmpfile)

    def test_check_illuminarun_complete(self):
        check_illuminarun_complete(
            'tests/data/ImageAnalysis_Netcopy_complete.txt', self.fp_tmpfile)
        with open(self.fp_tmpfile, 'r') as f:
            self.assertIn(
                "Raw Illumina data from run are complete.",
                '\n'.join(f.readlines()))

        with self.assertRaises(ValueError):
            check_illuminarun_complete(
                'tests/data/file_not_there.nowhere', self.fp_tmpfile)

        with self.assertRaises(ValueError):
            check_illuminarun_complete(
                'tests/data/ImageAnalysis_Netcopy_complete_fail.txt',
                self.fp_tmpfile)


if __name__ == '__main__':
    main()
