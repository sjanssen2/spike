from os.path import join, exists, dirname
from os import makedirs
import sys


def check_illuminarun_complete(fp_status_file, fp_output):
    """Test if illumina run data are complete by looking into a specific file.

    Parameters
    ----------
    fp_status_file : str
        Filepath of file from Illumina on which completion shall be tested.
        Currently <rundir>/ImageAnalysis_Netcopy_complete.txt.
    fp_output : str
        Filepath of reporting output file, if check is successful

    Raises
    ------
    ValueError if fp_status_file is not present or does not contain expected
    pattern.
    """
    PATTERN = ',Illumina RTA '
    run_complete = False

    if exists(fp_status_file):
        with open(fp_status_file, 'r') as f:
            for line in f.readlines():
                if PATTERN in line:
                    run_complete = True
                    break

    if run_complete is True:
        makedirs(dirname(fp_output), exist_ok=True)
        with open(fp_output, 'w') as f:
            f.write("Raw Illumina data from run are complete.\n")
    else:
        raise ValueError(
            ("According to file '%s' the raw Illumina data are not yet "
             "complete for this run.\nCheck file contents for "
             "pattern '%s'.\n") % (fp_status_file, PATTERN))
