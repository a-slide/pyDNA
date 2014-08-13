#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Standard library packages import
import csv

# Third party packages import
import pysam

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class PileUpMaker (object):
    """
    Class a generating a coverage graph from a dictionnary of coverage list indexed by sequence name
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__ (self, min_depth=1000, min_freq=0.02):
        """
        Create a PileUpMaker object
        @param min_depth Minimal depth to search for variations. Else the position will be ignored
        @param min_freq Mimimal frequency of a base at a given positionto be considered above the
        threshold. If bellow the base will be considered as undetected
        """
        # Creating object variables
        self.min_depth = min_depth
        self.min_freq = min_freq
        self.PileUp = ""
        
    def __repr__(self):
        msg = "PILEUP MAKER\n"
        msg+= "Minimal depth of positions : {}\n".format(self.min_depth)
        msg+= "Mimimal frequency of 'frequent' bases: {}\n".format(self.min_freq)
        return msg

    def __str__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)
        
    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def make (self, bam_path, outpath="./out", ref_name = "ref"):
        """
        """
        
        print ("\tCreate a pileUp file...")
        # Create a list and for results collecting
        out_list =[]
        
        # Open a handle on the bamfile and create an output file
        bamfile = pysam.Samfile(bam_path, "rb")
        # Open a PileUpColumn generator with an high limit of depth to avoid cuts
        # Analyse only if the sequencing depth is sufficient
        for PileUpCol in bamfile.pileup(max_depth=1000000):
            if PileUpCol.n > self.min_depth:

                # Fill a dictionary for positions of reads overlapping each positions.
                posDic = {'A':0,'C':0,'T':0,'G':0,'Del':0,'N':0}
                for read in PileUpCol.pileups:
                    if read.is_del:
                        posDic['Del']+=1
                    else:
                        try:
                            posDic[read.alignment.seq[read.qpos].upper()]+=1
                        except IndexError:
                            posDic['N']+=1

                # Define a threshold above which variations are not considered
                threshold = int(PileUpCol.n*self.min_freq)

                # If more than one frequent DNA base or indel was found at this position
                if sum([1 for base, count in posDic.items() if count >= threshold]) >= 2:
                    out_list.append([
                        ref_name,
                        bamfile.getrname(PileUpCol.tid),
                        PileUpCol.pos,
                        PileUpCol.n,
                        posDic['A'],
                        posDic['T'],
                        posDic['C'],
                        posDic['G'],
                        posDic['N'],
                        posDic['Del'],
                        round(posDic['A']/float(PileUpCol.n), 3),
                        round(posDic['T']/float(PileUpCol.n), 3),
                        round(posDic['C']/float(PileUpCol.n), 3),
                        round(posDic['G']/float(PileUpCol.n), 3),
                        round(posDic['N']/float(PileUpCol.n), 3),
                        round(posDic['Del']/float(PileUpCol.n), 3)])

        # Create a file to write out the list if results were found
        if out_list:
            self.PileUp =  "{}_{}_pileup.csv".format(outpath, ref_name)
            with open(self.PileUp, 'wb') as csvfile:
                writer = csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
                writer.writerow(["Ref","Seq","Pos","Total","CountA","CountT","CountC","CountG",
                    "CountN", "CountDel","FreqA","FreqT","FreqC","FreqG","FreqN","FreqDel"])
                for line in out_list:
                    writer.writerow(line)

        else:
            print("\t  No frequent variation found")
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class PileUpDecoy(object):
    """
    @class PileUpDecoy
    @brief Decoy class implementing no behaviour
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    def __init__ (self, *args, **kwargs):
        """
        Decoy init method
        """
        pass

    def __repr__(self):
        msg = "PILEUP DECOY\n"
        return msg

    def __str__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)
        
    def make (self, *args, **kwargs):
        """
        Decoy make method
        """
        print ("\tNo PileUp file to be generated")
