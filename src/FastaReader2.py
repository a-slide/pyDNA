#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages import
import gzip
from os import path, remove

# Third party packages import
from Bio import SeqIO

# Local Package import
from Utilities import file_extension, file_basename

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class FastaReader2(object):
    """
    @class  FastaReader2
    @brief
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        msg = self.__str__()
        msg += "List of references and sequences\n"
        for rname, seqlist in self.ref_seqname.items():
            msg += "{} : ".format (rname)
            msg += "{}\n\n".format ("  ".join(seqlist))
        if self.write_merge:
            msg += "Merged reference : {}\n".format (self.out)
        return msg

    def __str__(self):
        return "\n<Instance of {} from {}>\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, ref_list, write_merge=False, compressed_output=False, output="merged.fa"):
        """
        @param ref_list List of path to all ref files
        @param write_merge If True a merged fasta file will be generated
        @param compressed_output if True the output file will be gziped (much longer)
        @param output path to the output file
        """
        self.ref_list = ref_list
        self.write_merge = write_merge
        self.ref_seqname = {} # Will store refnames and associated seqnames

        if self.write_merge:
            self.compressed_output = compressed_output
            if self.compressed_output:
                self.out = path.abspath(output+".gz")
            else:
                self.out = path.abspath(output)
            self._read_and_merge()
        else:
            self._simple_read()

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _read_and_merge(self):
        """
        Read a list of fasta files, catches sequence names and write all sequences in the same
        output file
        """
        try:
            outfile = gzip.open(self.out,"wb") if self.compressed_output else open(self.out, "wb")

            # For all ref fasta file in the ref list
            for ref in self.ref_list:
                print ("Start parsing {}...".format(ref))

                # Open gzipped fasta or uncompressed fasta
                compressed_input = (ref[-2:].lower() == "gz")
                infile = gzip.open(ref,"rb") if compressed_input else open(ref,"rb")

                ref_name = file_basename (ref)
                self.ref_seqname[ref_name] = []

                # Iterate sequence with SeqIO parse
                for seq in SeqIO.parse(infile,"fasta"):
                    print seq
                    self._check_duplicate(seq.id, ref_name)
                    self.ref_seqname[ref_name].append(seq.id)
                    # And in every case write the line in the file
                    outfile.write(seq.format("fasta"))
                seq=None

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
                print ("Start parsing {}...".format(ref))

                # Open gzipped fasta or uncompressed fasta
                compressed_input = (ref[-2:].lower() == "gz")
                infile = gzip.open(ref,"rb") if compressed_input else open(ref,"rb")

                ref_name = file_basename (ref)
                self.ref_seqname[ref_name] = []

                # Iterate over sequence with SeqIO parse
                for seq in SeqIO.parse(infile, "fasta"):
                    self._check_duplicate(seq.id, ref_name)
                    self.ref_seqname[ref_name].append(seq.id)
                seq=None

                # Close the infile at the end of its parsing
                infile.close()

        # Catch all possible exception to close and infile
        except Exception as E:
            infile.close()
            raise Exception (E.message+"Can not create a list of sequence in the list of ref")


    def _check_duplicate (self, seqname, ref_name):
        """
        check if the entry is duplicated in the other list of self.ref_seqname
        """
        for rname, seqlist in self.ref_seqname.items():
            if seqname in seqlist:
                raise Exception ("{} from reference {} was already found in reference {}\n".format(
                    seqname, ref_name, rname))

