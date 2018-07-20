import sys
sys.path.append("..")

from unittest import TestCase, main
import tempfile
import os
from scripts.reports import *

class ChecksTests(TestCase):
    fp_tmpfile = None

    def setUp(self):
        _, self.fp_tmpfile = tempfile.mkstemp()

    def tearDown(self):
        os.remove(self.fp_tmpfile)

    def test_report_undertermined_filesizes(self):
        report_undertermined_filesizes('/home/jansses/spike2/data/aggregation/180614_SN737_0438_BCC7MCACXX.txt')


if __name__ == '__main__':
    main()
