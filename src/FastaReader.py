#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages import
import gzip
from tempfile import mkstemp
from os.path import getsize, isfile
from sys import exit as sys_exit
from os import remove

# Third party package import
from Bio import SeqIO

# Local Package import
from Utilities import file_extension, file_basename

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class FastaReader(object):
    """
    @class  FastaReader
    @brief
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        msg = self.__str__()
        for source in self.sources:
            msg += repr(source)
        return msg

    def __str__(self):
        return "\n<Instance of {} from {}>\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, source_paths, write_merge=False, output="merged.fa"):
        """
        @param source_paths List of path to all source files
        @param write_merge If True a merged fasta file containing both references will be
        generated in a temporary file
        @param
        """
        self.source_paths = source_paths
        self.sources = []

        if write_merge:
            self.output = output
            self._read_and_merge()
        else:
            self._simple_read()

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _read_and_merge(self):
        """
        """
        try:
            outfile = open(self.output, "w")

            # For all source fasta file in the source list
            for source in self.source_paths:
                print ("Start parsing {}...".format(source))
                # Open gzipped fasta or uncompressed fasta
                if file_extension(source)  in ["gz","GZ"]:
                    infile = gzip.open(source, "r")
                else:
                    infile = open(source, "r")

                # Parse the file
                seq_dict = {}
                for sequence in SeqIO.parse(infile, "fasta"):
                    # Fill the list and write the sequence in the merged reference
                    # Raise an error in case of redondancy
                    if sequence.id in seq_dict:
                        raise Exception ("Duplicated entry in {}".format(source))
                    else:
                        seq_dict[sequence.id] = len(sequence)
                        outfile.write(sequence.format("fasta"))
                #Creating a source object and adding it to the list
                self.sources.append(Source(source, file_basename(source), seq_dict))

            # Close and verify the output file
            outfile.close()
            # Verify if the list was filled if the file was generated
            assert self.sources, "No sequences were found in the source files"
            assert isfile(self.output), "No output file found"
            assert getsize(self.output) > 0, "The output file is empty"

        except Exception as E:
            print(E)
            print("Removing temporary file")
            try:
                if isfile (self.output):
                    remove (self.output)
            except Exception:
                pass
            exit

        return True


    def _simple_read(self):
        """
        Read fasta file sources and store names of sequences in lists
        """
        try:
             # For all source fasta file in the source list
            for source in self.source_paths:
                print ("Start parsing {}...".format(source))
                # Open gzipped fasta or uncompressed fasta
                if file_extension(source) == "gz" or file_extension(source) == "GZ":
                    infile = gzip.open(source, "r")
                else:
                    infile = open(source, "r")

                # Parse the file
                seq_dict = {}
                for sequence in SeqIO.parse(infile, "fasta"):
                    # Fill the list and write the sequence in the merged reference
                    # Raise an error in case of redondancy
                    if sequence.id in seq_dict:
                        raise Exception ("Duplicated entry in {}".format(source))
                    else:
                        seq_dict[sequence.id] = len(sequence)
                #Creating a source object and adding it to the list
                self.sources.append(Source(source, file_basename(source), seq_dict))

            # Verify if the list was filled if the file was generated
            assert self.sources, "No sequences were found in the source files"

        except Exception as E:
            print(E)
            exit

        return True

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Source (object):
    """
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    def __repr__(self):
        msg = "{} (@ {})\n".format(self.name, self.fasta_path)
        for name, size in self.seq_dict.items():
            msg += "\tSeq:{}\tLen:{}\n".format(name, size)
        return msg

    def __str__(self):
        return "\n<Instance of {} from {}>\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, fasta_path, name, seq_dict):
        """
        @param path Path of the fasta file containing the Source
        @param name Basename of the source
        @param seq_list Dictionnary of {name:size} of each sequence in the source
        """
        self.fasta_path = fasta_path
        self.name = name
        self.seq_dict = seq_dict

    def get_seq_names (self):
        return [name for name in self.seq_dict.keys()]

    def get_seq_sizes (self):
        return [size for size in self.seq_dict.values()]

    def get_seq (self):
        return [[name, size] for name, size in self.seq_dict.items()]
