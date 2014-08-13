#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Local Package import

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class BedGraphMaker (object):
    """
    Class a generating a coverage graph from a dictionnary of coverage list indexed by sequence name
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__ (self, min_depth=0, make_bed=False):
        """
        Create a BedGraphMaker object
        @param min_depth Minimal depth to report. Elsewhere it will be considered to be null
        @param make_bed Report all positions in the file even if bellow the minimal depth (= 0)
        """
        # Creating object variables
        self.min_depth = min_depth
        self.make_bed = make_bed
        self.bedgraph = ""
        self.bed = ""
        
    def __repr__(self):
        msg = "BEDGRAPH MAKER\n"
        msg+= "Minimal depth reported : {}\n".format(self.min_depth)
        msg+= "Create a Bed  report instead of a bedgraph : {}\n".format(self.bed)
        return msg

    def __str__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)
        
    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def make (self, cov_dict, outpath="./out", ref_name = "ref"):
        """
        @param cov_dict Dictionnary of coverage lists indexed by sequence name
        @param outpath Basename of the path where to output files
        @param ref_name Name of the reference genome containing the sequence listed in cov_dict
        """
        if self.make_bed:
            print ("\tCreate a bed File...")
            self.bed =  "{}_{}.bed".format(outpath, ref_name)
            with open (self.bed, "wb") as outfile:
                # Write bed header
                outfile.write ("track type=bed name={} color=0,0,0\n".format(ref_name))
                outfile.write(self._make_bed(cov_dict))
            
        else:
            print ("\tCreate a bedGraph File...")
            self.bedgraph =  "{}_{}.bedgraph".format(outpath, ref_name)
            with open (self.bedgraph, "wb") as outfile:
                # Write bedGraph header
                outfile.write ("track type=bedGraph name={} color=0,0,0\n".format(ref_name))
                outfile.write(self._make_bedgraph(cov_dict))

    def _make_bed (self, cov_dict):
        """
        Return a string for all positions in all sequences 
        """
        out = ""
        for seq_name, cov_list in cov_dict.items():
            for pos, depth in enumerate (cov_list):
                if depth < self.min_depth:
                    depth = 0
                out += ("{}\t{}\t{}\t{}\n".format(seq_name, pos, pos, depth))
        return out
                
    def _make_bedgraph (self, cov_dict):
        """
        Return a string for positions grouped by interval of the same depth in all sequences.
        only if it is above the threshold. 
        """
        out = ""
        for seq_name, cov_list in cov_dict.items():
            start = -1
            depth_prec = 0
            for pos, depth in enumerate (cov_list):
                if depth >= self.min_depth:
                    if depth != depth_prec:
                        if start != -1:
                            out += ("{}\t{}\t{}\t{}\n".format(seq_name, start, pos-1, depth_prec))
                        start = pos
                        depth_prec = depth

                else:
                    if start != -1:
                        out += ("{}\t{}\t{}\t{}\n".format(seq_name, start, pos-1, depth_prec))
                    start = -1
                    depth_prec = 0
        return out
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class BedGraphDecoy(object):
    """
    @class BedGraphDecoy
    @brief Decoy class implementing no behaviour
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    def __init__ (self, *args, **kwargs):
        """
        Decoy init method
        """
        pass

    def __repr__(self):
        msg = "BEDGRAPH DECOY\n"
        return msg

    def __str__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)
        
    def make (self, *args, **kwargs):
        """
        Decoy make method
        """
        print ("\tNo BedGraph file to be generated")
