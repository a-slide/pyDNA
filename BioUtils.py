# -*- coding: utf-8 -*-

"""
@package pyDNA
@brief  Collection of functions to manipulate DNA sequences
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2015
* <adrien.leger@gmail.com> <adrien.leger@inserm.fr> <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

# Standard library imports
from gzip import open as gopen
from unicodedata import normalize
from os import path
from glob import iglob
import csv

# Local imports
from FileUtils import is_gziped


def import_seq(filename, col_type="dict", seq_type="fasta"):
    """
    Import sequences from a fasta files in a list of biopython SeqRecord
    @param filename Valid path to a fasta file. Can contains several sequences and can be gzipped
    @param col_type Type of the collection where SeqReccord entries will be added ("list" or "dict").
    @param seq_type Type of the sequence file to parse (see Biopython seqIO for supported format)
    @return A list or a dictionnary containing all seqReccord objects from the fastq file
    @exception IOError  Raise if the path in invalid or unreadeable
    """
    # Specific Third party import
    from Bio import SeqIO

    # Try to open the file first gz compressed and uncompressed
    try:

        # Verify if the type of the input sequence is valid
        seq_type = seq_type.lower()
        allowed_seq = ["fasta", "genbank", "gb", "fastq-illumina", "fastq-solexa" , "fastq",
        "fastq-sanger", "embl ", "abi ", "seqxml", "sff", "uniprot-xml"]
        assert seq_type in allowed_seq , "The input file format have to be in the following list : "+ ", ".join(allowed_seq)

        # Verify if the type of the output collection is valid
        col_type = col_type.lower()
        allowed_types = ["dict", "list"]
        assert col_type  in allowed_types, "The output collection type have to be in the following list : "+ ", ".join(allowed_types)

        # Open gzipped or uncompressed file
        if is_gziped(filename):
            #print("\tUncompressing and extracting data")
            handle = gopen(filename, "r")
        else:
            #print("\tExtracting data")
            handle = open(filename, "r")

        # Create the collection
        if col_type == "list":
            seq_col = list(SeqIO.parse(handle, seq_type))
        else:
            seq_col = SeqIO.to_dict(SeqIO.parse(handle, seq_type))

        # Close file, verify if the collection is filled and returned it
        handle.close()
        assert seq_col, 'The collection contains no SeqRecord after file parsing. Exit'
        return seq_col

    except IOError as E:
        print(E)
        exit()

    except AssertionError as E:
        print (E)
        exit()

def count_seq (filename, seq_type="fasta"):
    """
    Count the number of sequences in a fastq or a fastq file
    @param filename Path to a valid readeable file
    @param file_type Should be either fastq or fastq. Default fasta
    """
    # Verify if the file is fasta or fastq type
    assert seq_type in ["fasta", "fastq"], "The file has to be either fastq or fasta format"

    # Open the file
    if is_gziped(filename):
        fp = gopen(filename, "rb")
    else:
        fp = open(filename, "rb")

    # line no counter
    nline = 0

    # FASTA Find a start line seq character ">" an increment the counter each time
    if seq_type ==  "fasta":
        for line in fp:
            if line[0] == ">":
                nline+=1
        fp.close()
        return nline

    # FASTQ No motif to find, but 4 lines correspond to 1 sequence
    else:
        for line in fp:
            nline+=1
        fp.close()
        return nline/4

def DNA_reverse_comp (sequence, AmbiguousBase=True):
    """
    Generate the reverese complementary sequence of a given DNA sequence
    @param sequence DNA sequence
    """
    if AmbiguousBase:
        compl = {'A':'T','T':'A','G':'C','C':'G','Y':'R','R':'Y','S':'S','W':'W','K':'M','M':'K',
               'B':'V','V':'B','D':'H','H':'D','N':'N'}
    else:
        compl = {'A':'T','T':'A','G':'C','C':'G','N':'N'}

    compl_sequence=""
    for base in sequence:
        try:
            compl_sequence += compl[base.upper()]
        except KeyError:
            compl_sequence += 'N'

    return compl_sequence[::-1]

def gb_to_bed(gb_file, track_description="User Supplied Track", features_type = []):
    """
    write a bed file containing annotations from a gb file
    @param gb_file path to a gb file containing a single molecule
    @param track_description Line of text descitpion for the bed track (60 char max)
    @param feature_type Restrict to the list of feature indicated
    """
    # Specific Third party import
    from Bio import SeqIO

    # open in and out files
    with open (gb_file, "rb") as gb:
        outf = file_basename(gb_file)+".bed"
        with open(outf, 'w') as bed:

            # parse record
            record = SeqIO.read( gb, "genbank")
            print ("{} features to be parsed\n".format(len(record.features)))

            # write bed header
            bed.write("""track name={} description="{}" visibility=2 colorByStrand="255,0,0 0,0,255"\n""".format(
                record.id,
                track_description))

            for feature in record.features:

                # if specific deatures are required
                if features_type and feature.type not in features_type:
                    continue

                start = feature.location.start.position
                stop = feature.location.end.position
                strand = '-' if feature.strand < 0 else '+'

                # try several alternative key for feature name else skip the feature
                if 'gene' in feature.qualifiers:
                    name = feature.qualifiers['gene'][0]
                elif 'locus_tag' in feature.qualifiers:
                    name = feature.qualifiers['locus_tag'][0]
                elif 'label' in feature.qualifiers:
                    name = feature.qualifiers['label'][0]
                else:
                    continue

                # Normalize name to asci
                name = str(name.decode('ascii', 'ignore'))
                bed.write("{0}\t{1}\t{2}\t{3}\t1000\t{4}\t\n".format(
                    record.id,
                    start,
                    stop,
                    name,
                    strand))
    return outf

def fetch_count_read (alignment_file, seq_name, start, end):
    """
    Count the number of read that are at least partly overlapping a specified chromosomic region
    @param alignment_file Path to a sam or a bam file
    @param seq_name Name of the sequence where read are to be aligned on
    @param start Start genomic coordinates of the area of alignment
    @param end End End genomic coordinates of the area of alignment
    """
    # Specific Third party import
    from pysam import AlignmentFile

    # Init a generator on the sam or bam file with pysam
    if alignment_file[-3:].lower() == "bam":
        al = AlignmentFile(alignment_file, "rb")

    elif alignment_file[-3:].lower() == "sam":
        al = AlignmentFile(alignment_file, "r")

    else:
        raise Exception("Wrong file format (sam or bam)")

    # Count read aligned at least partly on the specified region
    n = 0
    for i in al.fetch(seq_name, start, end):
        n += 1

    al.close()

    return n


def fetch_all_bam (bam_pattern, seq_name, coord_list, outname):
    """
    for all bam files matching the pattern, count the number of read overlapping list of coordinates
    and generate a file report
    @param bam_pattern Pattern to match in bam file to be included in the analysis
    @param seq_name Name of the sequence where read are to be aligned on
    @param coord_list list of coordinate start and end where to find overlapping reads
    @param outname Name of the output csv file repport
    """
    # Specific Third party import
    from pysam import index as pysamIndex

    # Find bam matching the patern and sorting the list alphabetically
    bam_list = list(iglob("*"+bam_pattern+"*"+".bam"))
    bam_list.sort()
    print ("Files analysed :")
    print (bam_list)

    # Generate a table header
    header = ["Coordinates"]
    for start, end in coord_list:
        header.append("{}:{}".format(start, end))

    # Generate a list to store hit found
    all_hits = []
    all_hits.append(header)

    # Fetch file for coordinates in coord list
    for bam in bam_list:
        print("\nAnalysing file : "+bam)

        if not path.isfile(bam+".bai"):
            print ("\tGenerate a Bam index")
            pysamIndex(bam)

        hits = [bam]
        for start, end in coord_list:
            hits.append(fetch_count_read(bam, seq_name, start, end))

        all_hits.append(hits)

    # Finally write a new table
    print ("\nWrite results in a csv table")
    with open(outname, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for i in all_hits:
            writer.writerow(i)

    print("Done")

    return
