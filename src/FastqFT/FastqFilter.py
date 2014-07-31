#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages import
import gzip
from time import time

# Third party package import
from Bio import SeqIO

# Local Package import
from Utilities import file_basename

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class FastqFilter(object):
    """
    @class  FastqFilter
    @brief Main class of the package
    Require the third party package Biopython
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        return self.__str__() + self.get_report()

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, quality_filter=None, adapter_trimmer=None, input_qual="fastq-sanger"):
        """
        """
        # Declare and initialize object variables
        self.total = 0

        # Store parameters in object variables
        self.qual = quality_filter
        self.adapt = adapter_trimmer
        self.input_qual = input_qual

    #~~~~~~~CLASS METHODS~~~~~~~#

    def filter(self, R1_path, R2_path):
        """
        """
        # Start a time counter
        start_time = time()

        try:
            print("Uncompressing and extracting data")
            # Open fastq and initialize a generator to iterate over files
            in_R1 = gzip.open(R1_path, "r")
            in_R2 = gzip.open(R2_path, "r")
            genR1 = SeqIO.parse(in_R1, self.input_qual)
            genR2 = SeqIO.parse(in_R2, self.input_qual)

            # Open output stream to write in fastq.gz files
            out_name_R1 = file_basename(R1_path) + "_trimmed.fastq.gz"
            out_name_R2 = file_basename(R2_path) + "_trimmed.fastq.gz"
            out_R1 = gzip.open(out_name_R1, "w")
            out_R2 = gzip.open(out_name_R2, "w")

            print("Parsing files and filtering sequences\n")
            # Parsing files and filtering sequences matching quality requirements
            while True:

                # Import a new pair of fastq until end of file
                seqR1 = next(genR1, None)
                seqR2 = next(genR2, None)
                if not seqR1 or not seqR2:
                    break
                self.total +=1
                if self.total%10000 == 0:
                    print ("\t{} sequences processed".format(self.total))

                # Quality filtering
                if self.qual:
                    seqR1 = self.qual.filter(seqR1)
                    seqR2 = self.qual.filter(seqR2)
                    if not seqR1 or not seqR2:
                        continue

                # Adapter trimming and size filtering
                if self.adapt:
                    seqR1 = self.adapt.trimmer(seqR1)
                    seqR2 = self.adapt.trimmer(seqR2)
                    if not seqR1 or not seqR2:
                        continue

                # If all test passed = write fastq to the files
                out_R1.write(seqR1.format("fastq-sanger"))
                out_R2.write(seqR2.format("fastq-sanger"))

            in_R1.close()
            in_R2.close()
            out_R1.close()
            out_R2.close()

        except IOError as E:
            print (E)
            exit (0)

        return (self.get_report(time()-start_time, R1_path, R2_path))


    def get_report (self, exec_time=None, R1_path=None, R2_path=None):
        """
        Return a report
        """
        report = "====== FASTQFILTER QUALITY CONTROL ======\n\n"
        if exec_time:
            report += "  Execution time : {} s\n".format(exec_time)
            report += "  R1 file : {}\n".format (R1_path)
            report += "  R2 file : {}\n".format (R2_path)
        report += "  Total sequences processed : {}\n".format(self.total)
        report += "  Input quality score : {}\n\n".format (self.input_qual)

        return report
