#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Local Package import
from Utilities import fill_between_graph

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class QualityFilter(object):
    """
    @class  QualityFilter
    @brief
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        msg = self.__str__()
        msg += "\tQuality Threshold : {}\n".format(self.min_qual)
        if self.run:
            msg += "\tTotal sequences : {}\n".format(self.total)
            msg += "\tFail quality filter : {}\n".format(self.qual_fail)
            msg += "\tPass quality filter : {}\n".format(self.qual_pass)
            msg += "\tMean quality : {}\n".format(sum(self.mean_qual)/float(len(self.mean_qual)))
            msg += "\tMinimal quality : {}\n".format(min(self.mean_qual))
            msg += "\tMaximal quality : {}\n".format(max(self.mean_qual))
        return msg

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__(self, min_qual):
        """
        """
        # Init object variables
        self.min_qual = min_qual
        self.total = 0
        self.qual_pass = 0
        self.qual_fail = 0
        self.mean_qual = []
        self.run = False

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def filter(self, record):
        """
        Compute mean quality score and compare to the minimal quality required
        """
        self.run = True
        # Compute the mean quality
        mean = sum(record.letter_annotations['phred_quality'])/len(record)
        # Add the value to the mean list
        self.mean_qual.append(mean)
        self.total += 1
        # Return the record if its quality is high enough
        if mean >= self.min_qual:
            self.qual_pass += 1
            return record
        else:
            self.qual_fail += 1
            return None

    #~~~~~~~ GETTERS ~~~~~~~#

    def get_mean_list (self):
        return self.mean_qual

    def get_mean_qual (self):
        if len(self.mean_qual) > 0 :
            return sum(self.mean_qual)/float(len(self.mean_qual))
        else:
            return None

    def get_min_qual (self):
        if len(self.mean_qual) > 0 :
            return min(self.mean_qual)
        else:
            return None

    def get_max_qual (self):
        if len(self.mean_qual) > 0 :
            return max(self.mean_qual)
        else:
            return None

    def get_tot_seq (self):
            return self.total

    def get_qual_fail (self):
            return self.qual_fail

    def get_qual_pass (self):
            return self.qual_pass
