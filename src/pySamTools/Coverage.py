#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Local Package import
from Utilities import fill_between_graph
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class CoverageMaker (object):
    """
    Class a generating a coverage graph from a dictionnary of coverage list indexed by sequence name
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__ (self, min_depth=0, make_bedgraph=True ,make_bed=False, make_covgraph=False):
        """
        Create a CoverageMaker object
        @param min_depth Minimal depth to report. Elsewhere it will be considered to be null
        @param make_bed Report all positions in the file even if bellow the minimal depth (= 0)
        """
        # Creating object variables
        self.min_depth = min_depth
        
        self.make_bedgraph = make_bedgraph
        self.bedgraph = ""
        
        self.make_bed = make_bed
        self.bed = ""

        self.make_covgraph = make_covgraph
        self.covgraph_list = []
        
    def __repr__(self):
        msg = "\tCOVERAGE MAKER\n"
        
        if not self.make_bed and not self.make_bedgraph and not self.make_covgraph:
            msg += "\t\tNo output requested\n"
            return msg
            
        msg+= "\t\tMinimal depth : {}\n".format(self.min_depth)
        msg+= "\t\tOutput requested :"
        if self.make_bedgraph:
            msg+= "\tBedGraph"
        if self.make_bed:
            msg+= "\tBed"
        if self.make_covgraph:
            msg+= "\tCovGraph"
        msg+= "\n"
        return msg

    def __str__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)
        
    def get(self, key):
        return self.__dict__[key]

    def set(self, key, value):
        self.__dict__[key] = value
        
    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def make (self, cov_dict, outpath="./out", ref_name = "ref"):
        """
        @param cov_dict Dictionnary of coverage lists indexed by sequence name
        @param outpath Basename of the path where to output files
        @param ref_name Name of the reference genome containing the sequence listed in cov_dict
        """
        
        if not self.make_bed and not self.make_bedgraph and not self.make_covgraph:
            return
        
        if self.make_bedgraph:
            print ("\tCreate a bedGraph File...")
            self.bedgraph =  "{}_{}.bedgraph".format(outpath, ref_name)
            with open (self.bedgraph, "wb") as outfile:
                # Write bedGraph header
                outfile.write ("track type=bedGraph name={} color=0,0,0\n".format(ref_name))
                outfile.write(self._make_bedgraph(cov_dict))
        
        if self.make_bed:
            print ("\tCreate a bed File...")
            self.bed =  "{}_{}.bed".format(outpath, ref_name)
            with open (self.bed, "wb") as outfile:
                # Write bed header
                outfile.write ("track type=bed name={} color=0,0,0\n".format(ref_name))
                outfile.write(self._make_bed(cov_dict))
        
        if self.make_covgraph:
            print ("\tCreate coverage graphics...")
            self._make_covgraph(cov_dict, outpath="./out", ref_name = "ref")
            
            
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
    
    
    def _make_covgraph (self, cov_dict, outpath="./out", ref_name = "ref"):
        """
        @param cov_dict Dictionnary of coverage lists indexed by sequence name
        @param outpath Basename of the path where to output files
        @param ref_name Name of the reference genome containing the sequence listed in cov_dict
        """
        for seq_name, cov_list in cov_dict.items():
            
            # Filter coverage if threshold is higher than 0
            if self.min_depth > 0:
                for i in cov_list:
                    if i < self.min_depth:
                        i = 0
            
            # Create the graph with Utilities.fill_between_graph
            try:
                fill_between_graph (
                    X = [i+1 for i in range (len(cov_list))],
                    Y = cov_list,
                    basename = "{}_{}_{}".format(outpath,ref_name,seq_name),
                    img_type = "svg",
                    title = ("Coverage of reads over {} from {}".format (seq_name, ref_name)),
                    xlabel = 'Position',
                    ylabel = 'Count',
                    xsize = 50,
                    ysize = 10,
                    dpi = 150)
                    
                self.covgraph_list.append("{}_{}_{}.svg".format(outpath,ref_name,seq_name))
            
            except ImportError as E:
                print(E)
                print("Cannot create the required CovGraph file. Skip to the next step")
    
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class CoverageDecoy(object):
    """
    Decoy class implementing no behaviour
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    def __init__ (self, *args, **kwargs):
        """
        Decoy init method
        """
        pass

    def __repr__(self):
        msg = "COVERAGE DECOY\n"
        return msg

    def __str__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)
        
    def make (self, *args, **kwargs):
        """
        Decoy make method
        """
        print ("\tNo file to be generated")
