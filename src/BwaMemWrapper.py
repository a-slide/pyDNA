"""
@package BwaMemWrapper
@brief **Wrapper for BWA mem**
Please see the BWA user manual document for further details
[MANUAL](http://bio-bwa.sourceforge.net/bwa.shtml)
* To use the wrapper, it is first needed to generate a bwa index either by using NewIndex if no
index was already created or ExistingIndex if an index is already available.
* A instance of BwaMem can then be created by passing the index object as an argument.
* The BwaMem method align can finally be used as many times as desired with different queries.
It will returned each time a path to a sam file

@copyright [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages import
from os import path, remove, rmdir
from multiprocessing import cpu_count
from time import time

# Local library packages
from Utilities import run_command, file_basename, make_cmd_str

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class BwaMem(object):
    """
    @class  BwaMem
    @brief  Perform a lignement of a query fastq file against a bwa index. Results are written in
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

    def __init__ (self, Index, align_opt=None, bwa_mem = "bwa mem"):
        """
        Initialize the object and index the reference genome if necessary
        @param Index Bwa index object NewIndex or ExistingIndex
        @param align_opt Bwa mem dictionnary of option arguments such as "-t 5". The option flag
        have to be the key (without "-") and the the option value in the dictionnary value. If no
        value is requested after the option flag "None" had to be asigned to the value field.
        @param bwa_mem Path ot the bwa executable. If bwa if already added to your system path
        do not change the default value
        """
        # Creating object variables
        self.aligner = bwa_mem
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class NewIndex(object):
    """
    @class NewIndex
    @brief Wrapper for bwa index. Create a reference index from a fasta file
    BWA 0.7.5+ needs to be install and eventually added to the path.
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        msg = self.__str__()
        msg += "bwa index path : {}\n".format(self.indexer)
        msg += "Blastn database path : {}\n".format(self.index_path)
        msg += "Options :\n"
        for i, j in self.index_opt.items():
            msg += "\tFlag : {}\tValue : {}\n".format(i,j)
        return msg

    def __str__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, ref_path, index_path="./out.idx", index_opt=None, bwa_index = "bwa index"):
        """
        Initialize the object and index the reference genome if necessary
        @param ref_path Path of the fasta file containing the reference sequence
        @param index_path Outname for the bwa index files basename.
        @param index_opt bwa index dictionnary of option arguments such as "-t 5". The option flag
        have to be the key (without "-") and the the option value in the dictionnary value. If no
        value is requested after the option flag "None" had to be asigned to the value field.
        @param bwa index Path ot the bwa executable. If bwa if already added to your system path
        do not change the default value
        """
        # Creating object variables
        self.indexer = bwa_index
        self.index_path = index_path
        self.ref_path = ref_path

        # init an option dict and attribute defaut options
        self.index_opt = index_opt if index_opt else {}
        if "a" not in self.index_opt:
            self.index_opt["a"] = "bwtsw"
        if "p" not in self.index_opt:
            self.index_opt["p"] = index_path

        # Make a db with makeblastdb if not possible = remove the files
        try:
            self._make_index()

        except Exception as E:
            self._remove_index_files()
            raise Exception (E.message+"Impossible to generate a valid index from the reference sequence")


    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _make_index(self):
        """
        Create a bwa index from ref_path using bwa index
        """
        # Build the command line thanks to make_cmd_str
        opt_list = [self.ref_path]
        cmd = make_cmd_str(self.indexer, self.index_opt, opt_list)
        print "Creating a BWA index for {}".format(file_basename(self.ref_path))


        # Run the command line without stdin and asking both stdout and stderr
        start_time = time()
        stdout, stderr = run_command(cmd, None, True, True)

        # Verify the output
        if not stdout:
            raise Exception ("Error, no data received from standard output\n"+stderr)

        print (stdout)
        print ("Index created in {}s".format(round(time()-start_time, 3)))

    def _remove_db_files(self):
        """
        Remove index files in case of exception during indexing
        """
        print("Remove index files")

        # Removing Index file
        for ext in ["amb", "ann", "bwt", "pac", "sa"]:
            f = "{}.{}".format(self.index_path, ext)
            if path.isfile (f):
                remove (f)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class ExistingIndex(object):
    """
    @class ExistingIndex
    @brief Import an existing bwa index + verify the existence of files
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    def __repr__(self):
        msg = self.__str__()
        msg += "Bwa Index Path : {}\n".format(self.index_path)
        return msg

    def __str__(self):
        return "\n<Instance of {} from {}>\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, index_path):
        """
        @param index_path The index path is the name of any of the index files up to but not
        including the final "amb", "ann", "bwt", "pac" and "sa"
        """
        # Creating object variables
        self.index_path = index_path

        print ("Checking index files")
        # Checking if all index files needed by bwa are
        for ext in ["amb", "ann", "bwt", "pac", "sa"]:
            f = "{}.{}".format(self.index_path, ext)

            if not path.isfile (f):
                raise Exception ("Invalid database : {} does not exist".format(f))
            if path.getsize(f) == 0:
                raise Exception ("Invalid database : {} is empty".format(f))

        print ("All index files are valid")
