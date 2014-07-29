"""
@package RefMasker
@brief Find DNA sequences homologies between a list of query and a subject and to write a modified
version of the subject sequence.
* First a blast database is created from a reference fasta file if it is not provided by the user.
* Then, query sequences fasta files are blasted against the newly created subject database.
* Finally, if blast hits are found (homologies) the subject genome is imported, hit locations are
masked with "Ns" and a new masked reference is written
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
from os import remove, path
import gzip

# Third party package import
from Bio import SeqIO

# Local library packages import
from Utilities import import_seq, file_basename, mkdir
from Blast import Blastn

#~~~~~~~MAIN METHODS~~~~~~~#

def mask (  subject_fasta,
            hit_list,
            ref_outdir="./references/",
            ref_outname="masked_ref.fa",
            compress_ouput=True ):
    """
    @return If the sequence was edited the path of the edited reference is indicated, else the
    unmodified sequence.
    """

    # Test if object the first object of hit_list have the require s_id, s_start and s_end fields
    try:
        a = hit_list[0].s_id
        a = hit_list[0].s_start
        a = hit_list[0].s_start

    except (IndexError, AttributeError) as E:
        print (E)
        print ("The hit_list does not contain suitable hit objects")
        print ("The subject fasta file will not be edited")
        return subject_fasta

    # Initialize output folder
    mkdir(ref_outdir)

    # Initialize input fasta file
    if subject_fasta[-2:].lower() == "gz":
        in_handle = gzip.open(subject_fasta, "r")
    else:
        in_handle = open(subject_fasta, "r")

    # Initialize output fasta file
    if compress_ouput:
        ref_path = path.join (ref_outdir, ref_outname+".gz")
        out_handle = gzip.open(ref_path, 'w')
    else:
        ref_path = path.join (ref_outdir, ref_outname)
        out_handle = open(ref_path, 'w')

    # Generate a list of ref that will need to be modified
    id_list = {hit.s_id:0 for hit in hit_list}.keys()

    # Iterate over record in the subject fasta file
    print ("\nWriting new version of {} in which homologies with the hit list are masked".format(
        file_basename (subject_fasta)))
    # Iterate over record in the subject fasta file
    for record in SeqIO.parse(in_handle, "fasta"):

        # Check if the record is in the list of record to modify
        if record.id in id_list:

            print ("Hit found in {}. Editing the sequence".format(record.id))
            # Casting Seq type to MutableSeq Type to allow string editing
            record.seq = record.seq.tomutable()

            # For each hit in the list of hit found
            for hit in hit_list:
                if record.id == hit.s_id:

                    # For all position between start and end coordinates modify the base by N
                    for position in range (hit.s_start, hit.s_end):
                        record.seq[position]= 'N'
        else:
            print ("No hit found in {}".format(record.id))

        # Finally write the sequence modified or not
        out_handle.write(record.format("fasta"))

    # close files and return the masked ref path
    in_handle.close()
    out_handle.close()
    return ref_path
