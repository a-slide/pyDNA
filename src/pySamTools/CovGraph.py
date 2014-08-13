#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Local Package import
from Utilities import fill_between_graph

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class CovGraphMaker (object):
    """
    Class a generating a coverage graph from a dictionnary of coverage list indexed by sequence name
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__ (self, min_depth=0):
        """
        Create a CovGraphMaker object
        @param min_depth Minimal depth to report. Elsewhere it will be considered to be null
        """
        # Creating object variables
        self.min_depth = min_depth
        self.graph_list =[]
        
    def __repr__(self):
        msg = "COVGRAPH MAKER\n"
        msg+= "Minimal depth reported : {}\n".format(self.min_depth)
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
        
        print ("\tCreate coverage graphics...")
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
                    
                self.graph_list.append("{}_{}_{}.svg".format(outpath,ref_name,seq_name))
            
            except ImportError as E:
                print(E)
                print("Cannot create the required CovGraph file. Skip to the next step")
            
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class CovGraphDecoy(object):
    """
    @class CovGraphDecoy
    @brief Decoy class implementing no behaviour
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    def __init__ (self, *args, **kwargs):
        """
        Decoy init method
        """
        pass

    def __repr__(self):
        msg = "COVGRAPH DECOY\n"
        return msg

    def __str__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)
        
    def make (self, *args, **kwargs):
        """
        Decoy make method
        """
        print ("\tNo CovGraph file to be generated")
