#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages import
import gzip
from os import path, remove
from sys import stdout

# Third party packages import
from Bio import SeqIO

# Local Package import
from Utilities import file_extension, file_name, file_basename

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class FastaReader2(object):
    """
    @class  FastaReader2
    @brief  Parse a list of fasta files and extract sequence names from fasta header.
    Additionnaly a merged reference fasta file containing all references from the list can be
    written on disk ( gz compressed or not). DUPLICATED SEQUENCE NAMES ARE NOT ALLOWED and will
    raise an exception.
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        msg = self.__str__()
        msg += "List of references and sequences\n"
        for rname, seqlist in self.ref_seqname.items():
            msg += "{} : ".format (rname)
            msg += "{}\n".format ("  ".join(seqlist))
        msg += "\n"
        if self.write_merge:
            msg += "Merged reference : {}\n".format (self.out)
        if self.write_summary:
            msg += "Merged summary : {}\n".format (self.outsum)

        return msg

    def __str__(self):
        return "\n<Instance of {} from {}>\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self,
        ref_list,
        write_merge=False,
        compress_output=False,
        output="./merge_ref.fa",
        write_summary=False):
        """
        @param ref_list List of path to all ref files
        @param write_merge If True a merged fasta file will be generated
        @param compress_output If True the output file will be gziped (much longer)
        @param output Path to the output file
        @param write_summary if true a simple text file representing self.ref_seqname will be created
        """
        # Store object variables
        self.ref_list = ref_list
        self.write_merge = write_merge
        self.write_summary = write_summary
        self.ref_seqname = {} # Will store refnames and associated seqnames

        # Store additional variables and call the suitabe parser according to write_merge Flag
        if self.write_merge:
            self.compress_output = compress_output
            self.out = path.abspath(output+".gz") if compress_output else path.abspath(output)
            self._read_and_merge()
        else:
            self._simple_read()

        if self.write_summary:
            self.outsum = path.abspath(output+".txt")
            self._write_summary()

    #~~~~~~~GETERS~~~~~~~#

    def getSeqFromRef (self, refname):
        return self.ref_seqname[refname]

    def getAllSeq (self):
        all_seqlist=[]
        for seqlist in self.ref_seqname.values():
            all_seqlist.extend(seqlist)
        return all_seqlist

    def seqOrigin (self, seqname):
        for rname, seqlist in self.ref_seqname.items():
            if seqname in seqlist:
                return rname
        return None

    def getMergeRef (self):
        return self.out

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _read_and_merge(self):
        """
        Read a list of fasta files, catches sequence names and write all sequences in the same
        output file. additionnaly write a simple descriptive text file containing the names of
        the
        """
        try:
            outfile = gzip.open(self.out,"wb") if self.compress_output else open(self.out, "wb")

            # For all ref fasta file in the ref list
            for ref in self.ref_list:
                stdout.write("Start parsing {} ".format(file_name(ref)))
                stdout.flush()

                # Open gzipped fasta or uncompressed fasta
                infile = gzip.open(ref,"rb") if (ref[-2:].lower() == "gz") else open(ref,"rb")

                ref_name = file_basename (ref)
                self.ref_seqname[ref_name] = []

                # Iterate sequence with SeqIO parse
                for seq in SeqIO.parse(infile,"fasta"):
                    self._check_dup(seq.id, ref_name)
                    self.ref_seqname[ref_name].append(seq.id)
                    # And in every case write the line in the file
                    outfile.write(seq.format("fasta"))
                    stdout.write("*")
                    stdout.flush()
                seq=None
                print ("")

                # Close the infile at the end of its parsing
                infile.close()

            # Close files
            outfile.close()

        # Catch all possible exception to close and remove files properly
        except Exception as E:
            infile.close()
            outfile.close()
            remove(self.out)
            raise Exception (E.message+"Can not create a merged reference")


    def _simple_read (self):
        """
        Read a list of fasta files, catches sequence names and write all sequences in the same
        output file
        """
        try:
            # For all ref fasta file in the ref list
            for ref in self.ref_list:
                stdout.write("Start parsing {} ".format(file_name(ref)))
                stdout.flush()

                # Open gzipped fasta or uncompressed fasta
                infile = gzip.open(ref,"rb") if (ref[-2:].lower() == "gz") else open(ref,"rb")

                ref_name = file_basename (ref)
                self.ref_seqname[ref_name] = []

                # Iterate over sequence with SeqIO parse
                for seq in SeqIO.parse(infile, "fasta"):
                    self._check_dup(seq.id, ref_name)
                    self.ref_seqname[ref_name].append(seq.id)
                    stdout.write("*")
                    stdout.flush()
                seq=None
                print ("")

                # Close the infile at the end of its parsing
                infile.close()

        # Catch all possible exception to close and infile
        except Exception as E:
            infile.close()
            raise Exception (E.message+"Can not create a list of sequence in the list of ref")


    def _check_dup (self, seqname, ref_name):
        """
        check if the entry is duplicated in the other list of self.ref_seqname
        """
        for rname, seqlist in self.ref_seqname.items():
            if seqname in seqlist:
                raise Exception ("{} from reference {} was already found in reference {}\n".format(
                    seqname, ref_name, rname))

    def _write_summary (self):
        """
        write a simple text file representing self.ref_seqname
        """
        with open(self.outsum, "wb") as outfile:
            for rname, seqlist in self.ref_seqname.items():
                outfile.write(rname+"\n")
                for seq in seqlist:
                    outfile.write(" "+seq+"\n")
