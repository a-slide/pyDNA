#~~~~~~~GLOBAL IMPORTS~~~~~~~#

# Local Package import
from Utilities import import_seq

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class AdapterTrimmer(object):
    """
    @class  AdapterTrimmer
    @brief
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        msg = "ADAPTER TRIMMER\n"
        msg += "  List of adapters imported for trimming"

        for a in self.adapter_list:
            msg += "\n  Name : {}".format(a.id)
            msg += "  Sequence : {}".format(a.seq)
            msg += "  Min score : {}".format(a.annotations["min_score"])
            msg += "  Min len : {}".format(a.annotations["min_len"])
            if self.run:
                msg += "  Trimmed : {}".format(a.annotations["count"])

        if self.run:
            msg += "\n  Sequences untrimmed : {}\n".format(self.seq_untrimmed)
            msg += "  Sequences trimmed : {}\n".format(self.seq_trimmed)
            msg += "  DNA base trimmed : {}\n".format(self.base_trimmed)
            msg += "  Fail len filtering: {}\n".format(self.len_fail)
            msg += "  Pass len filtering : {}\n".format(self.len_pass)
            msg += "  Total pass : {}\n\n".format(self.len_pass+self.seq_untrimmed)
        return msg

    def __str__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, aligner, adapter_path, min_read_len=0.6, min_match_len=0.8, min_match_score=1.4):
        """
        @param aligner Wrapper object of pairwise alignement
        @param adapter_path Path of a multi fasta file containing all the adapters sequences
        @param min_read_len Fraction of read lenth = minimal size of fragment after trimming
        @param min_match_len Minimal fraction of adapter len that needs to be aligned on the target
        @param min_match_score Minimal score per base for the alignment of adapter and read
        """
        # Init object variables
        self.min_read_len = min_read_len
        self.min_match_len = min_match_len
        self.min_match_score = min_match_score
        # wrapper of pairwise alignement
        self.aligner = aligner

        # Import a list of adapters and add the reverse complements of adapters to the list
        self.adapter_list = import_seq(adapter_path, "list", "fasta")
        adapter_rc = []
        for i in self.adapter_list:
            rc = i.reverse_complement()
            rc.id = i.id + "_RC"
            adapter_rc.append(rc)
        self.adapter_list.extend(adapter_rc)

        # Initialize entries in each adapter annotations dictionnaries
        # Count = counter of match
        # Min_score = minimal score of match normalized by the size of the adapter
        # Min_len = minimal size of the adapter to be aligned on the read
        for a in self.adapter_list:
            a.annotations["count"] = 0
            a.annotations["min_score"] = int(self.min_match_score * len(a))
            a.annotations["min_len"] = int(self.min_match_len * len(a))

        # Initialize generic counters
        self.seq_untrimmed = 0
        self.seq_trimmed = 0
        self.base_trimmed = 0
        self.len_pass = 0
        self.len_fail = 0
        self.run = False

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def trimmer (self, record):
        """
        @param record A seq record object containing the subject reference sequence to be trimmed
        """
        self.run = True
        match_list = []
        len_rec = len(record)

        # Set a new reference sequence into the aligner
        self.aligner.set_ref(str(record.seq))

        # Iterate over the adaters in adapter_list
        for a in self.adapter_list:

            # Find match of the adapter along the current read
            match = self.aligner.align(str(a.seq), a.annotations["min_score"], a.annotations["min_len"])

            # if a match was found = increment the counter and append the match to the match list
            if match:
                #print ("Adapter found : {}\tScore : {}\tCigar : {}\tMatchLen : {}".format(
                #a.id, match.score, match.cigar_string, match.query_end-match.query_begin))
                a.annotations["count"] += 1
                match_list.append(match)

        # In case no match were found, the sequence doesn't need to be modify
        if not match_list:
            self.seq_untrimmed += 1
            return record

        # Else find the longer interval without adaptor matches
        start, end = self._longer_interval (match_list, len_rec)

        # Update counters
        self.seq_trimmed += 1
        self.base_trimmed += (len_rec-(end-start))

        # Return a slice of the reccord corresponding to the longer interval
        if end-start >= int(self.min_read_len*len_rec):
            self.len_pass +=1
            return record[start:end]
        # Or None if smaller than min_size
        else:
            self.len_fail +=1
            return None

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _longer_interval(self, match_list, len_seq):
        """
        Find the first larger interval that do not overlapp any match in match list
        This strategy allow to use an unsorted list of match
        """
        # Initialize a list of boolean to False of the same size as the read
        coverage = [False for i in range(len_seq)]

        # Flag positions overlapped by a read by changing the boolean to True
        for m in match_list:
            for i in range (m.ref_begin, m.ref_end):
                coverage[i] = True

        # Read through the list to find the longer inteval between True flags
        start_max = end_max = inter_max = start = inter = 0
        for i in range(len_seq):
            if coverage[i]:
                start = i+1
                inter = 0
            else:
                inter += 1
                if inter > inter_max:
                    inter_max = inter
                    start_max = start
                    end_max = i

        #print ("Longer interval = {} [{}:{}]".format(inter_max, start_max+1, end_max-1))
        return start_max, end_max

    #~~~~~~~ GETTERS ~~~~~~~#

    def get_seq_untrimmed (self):
        return self.seq_untrimmed

    def get_seq_trimmed (self):
        return self.seq_trimmed

    def get_base_trimmed (self):
        return self.base_trimmed

    def get_len_pass (self):
        return self.len_pass

    def get_len_fail (self):
        return self.len_fail


