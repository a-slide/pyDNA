#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Standard library packages import
from os import remove, path
from tempfile import mkstemp
import gzip

# Local library packages import
from Utilities import mkdir, import_seq, file_basename, file_name, file_extension, fgunzip
from BlastnWrapper import Aligner
from MakeblastdbWrapper import NewDB, ExistingDB


#~~~~~~~MAIN METHODS~~~~~~~#

def align  (query_list,
            subject_db = None,
            subject_fasta = None,
            evalue=0.1,
            aligner = "blastn",
            align_opt=None,
            db_maker = "makeblastdb",
            db_opt=None,
            db_outdir = "./blast_db/",
            db_outname = "out"):

    """
    Main function of RefMasker that integrate database creation, blast and homology masking
    * Instantiate Blast database and blastn object
    * Perform iterative blasts of query sequences against the subject database and create a list of
    hits.
    @param query_list List of paths indicating fasta files containing query sequences (can be
    gzipped). Fasta can contains multiple sequences.
    @param subject_db Basename of file from a blast database created by "makeblastdb" if available
    @param subject_fasta Reference fasta file. Required if no ref_index is given (can be gzipped)
    @param evalue  Cutoff used in blast to select valid hits
    @param aligner Path ot the blastn executable. Not required if blast+ if added to your path
    @param blastn_opt Dict of options for blastn (see Utilities.make_cmd_str)
    @param db_maker Path ot the makeblastdb executable. Not required if blast+ if added to your path
    @param db_opt Dict of options for makeblastdb (see Utilities.make_cmd_str)
    @param db_outdir Directory where to store the database files
    @param db_outname Basename of the database files
    @return A list of BlastHit objects
    """
    # Try to import an existing database
    try:
        if not subject_db:
            raise Exception("No Blast database was provided")

        print("Existing database provided")
        db = ExistingDB(subject_db)

    # If no DB or if an error occured during validation of the existing DB = create a new db
    except Exception as E:
        print (E)

        # Verify the presence of the reference fasta file
        if not subject_fasta or not path.isfile (subject_fasta):
            raise Exception("Invalid or no fasta file provided. Cannot create a database")

        print ("Generate a database...")
        mkdir(db_outdir)
        db_path = path.join (db_outdir, db_outname)

        # If the fastq file is compressed = extract the file in a tempory file
        if file_extension(subject_fasta).lower() == "gz":
            tmp_handle, tmp_path = mkstemp(suffix='.fa')
            fgunzip (in_path=subject_fasta, out_path=tmp_path)
            db = NewDB(tmp_path, db_path, db_opt, db_maker)
            remove(tmp_path)
        else:
            db = NewDB(subject_fasta, db_path, db_opt, db_maker)

    # Initialise a Blastn object
    blast = Aligner(db, align_opt, aligner)

    # Generate a list of hit containing hits of all sequence in query list in subject
    hit_list = []
    # Extend the list of hits for each query in a bigger list.
    for query in query_list:
        hit_list.extend(blast.align(query, evalue))

    return hit_list
