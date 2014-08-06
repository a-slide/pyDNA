#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages import
from os import path, remove
from time import time

# Local library packages
from Utilities import run_command, file_basename, make_cmd_str

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
        msg = "BWA INDEX WRAPPER (NEW INDEX)\n"
        msg += "bwa index path : {}\n".format(self.indexer)
        msg += "Blastn database path : {}\n".format(self.index_path)
        msg += "Options : {}\n".format(self.index_opt)
        return msg

    def __str__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, ref_path, index_path="./out.idx", index_opt="", bwa_index = "bwa index"):
        """
        Initialize the object and index the reference genome if necessary
        @param ref_path Path of the fasta file containing the reference sequence
        @param index_path Outname for the bwa index files basename.
        @param index_opt bwa index dictionnary of option arguments such as "-t 5". The option flag
        have to be the key (without "-") and the the option value in the dictionnary value. If no
        value is requested after the option flag "None" had to be asigned to the value field.
        @param bwa index Path ot the bwa index executable. Not required if bwa if added to your path
        """
        # Creating object variables
        self.indexer = bwa_index
        self.index_path = index_path
        self.ref_path = ref_path

        # init an option string and attribute defaut options
        self.index_opt = "{} -a {} -p {}".format (index_opt, "bwtsw", self.index_path)
        
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
        # Build the command line
        cmd = "{} {} {}".format(self.indexer, self.index_opt, self.ref_path)
        
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
        msg = "BWA INDEX WRAPPER (EXISTING INDEX)\n"
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
