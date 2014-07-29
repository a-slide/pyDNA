#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages import
from os import path, remove, rmdir
from multiprocessing import cpu_count
from time import time

# Local library packages
from Utilities import run_command, file_basename, make_cmd_str

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Aligner(object):
    """
    @class  Aligner
    @brief  Perform a alignement of a query fastq file against a bwa index. Results are written in
    a sam file whose path is returned at the end of the alignment
    BWA 0.7.5+ needs to be install and eventually added to the path.
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        msg = self.__str__()
        msg += "Bwa mem path : {}\n".format(self.aligner)
        msg += "Options :\n"
        for i, j in self.align_opt.items():
            msg += "\tFlag : {}\tValue : {}\n".format(i,j)
        msg += "BwaIndex : {}\n".format(repr(self.Index))
        return msg

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, Index, align_opt=None, aligner = "bwa mem"):
        """
        Initialize the object and index the reference genome if necessary
        @param Index Bwa index object NewIndex or ExistingIndex
        @param align_opt Bwa mem dictionnary of option arguments such as "-t 5". The option flag
        have to be the key (without "-") and the the option value in the dictionnary value. If no
        value is requested after the option flag "None" had to be asigned to the value field.
        @param bwa_mem Path ot the bwa mem executable. Not required if bwa if added to your path
        """
        # Creating object variables
        self.aligner = aligner
        self.Index = Index
        self.align_opt = align_opt if align_opt else {}

        # By default the option t (number of thread to use) is setted to the max number of
        # available threads
        if "t" not in self.align_opt:
            self.align_opt["t"] = cpu_count()

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def align(self, R1_path, R2_path=None, out_path="./out.sam"):
        """
        Align query fastq against a subject database and return a list of BlastHit object
        @param R1_path Path to the file containing fastq sequences
        @param R2_path Facultative path to the file containing paired fastq sequence
        @param out_path Path to the output sam file
        @return A list of BlastHit objects if at least one hit was found
        @exception (SystemError,OSerror) May be returned by run_command in case of invalid command line.
        """
        # Build the command line with make_cmd_str
        # Pair end mode
        if R2_path:
            opt_list = [self.Index.index_path, R1_path, R2_path, "> "+out_path]
        # Single end mode
        else:
            opt_list = [self.Index.index_path, R1_path, "> "+out_path]
        cmd = make_cmd_str(self.aligner, self.align_opt, opt_list)

        # Execute bwa mem (Can raise a SystemError) and verify if stdout is not None
        print ("Align against {} index with bwa mem".format(file_basename (self.Index.index_path)))
        stderr = run_command(cmd, stdin=None, ret_stderr=True, ret_stdout=False)

        # In bwa stderr return a report of alignment
        print (stderr)
        return out_path
