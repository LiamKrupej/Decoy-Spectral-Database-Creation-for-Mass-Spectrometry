import sys, os
from create_spectral_decoys.custom_filtering import find_chem_string, are_peaks_similar, are_spectrums_same
from matplotlib import pyplot as plt
from analysis import Analysis
import numpy

# class implements Analysis class

class DecoyAnalysis(Analysis):
    
        # elements include a list of sorted scores, total false positives and total true hits
        # FDR list, qvalues, and a label. PIT set to 0 as not used in method implementation.
        
    def __init__(self, sorted_scores, PIT, label):
        self.sorted_scores = sorted_scores
        self.PIT = PIT
        self.libscores = []
        self.decoyscores = []
        self.fdrlist = []
        self.qvalues = []
        self.label = label
        
        for s in sorted_scores:
            if(s.type == "library"):
                self.libscores.append(s.score)
            else:
                self.decoyscores.append(s.score)
                
    # implementation for reversing a sequence of objects
    
    def enumerate_reversed(self, seq):
        n = len(seq)
        for obj in reversed(seq):
            n -= 1
            yield n, obj        
    
    # creates a list of fdr values using the sorted scores 
    
    
    def create_fdr_list(self): 
        if self.fdrlist:
            print("FDRs are already created for this set")
            return False
        
        else:
            
            print("PIT is ", self.PIT)
            for score in self.sorted_scores:
                if(score.type == "library"):
                    list_to_consider = [x for x in self.sorted_scores if x.score >= score.score]
                    decoy_hits_in_list = 0
                    library_hits_in_list = 0
                    for s in list_to_consider:     
                        if(s.type == "library"):
                            library_hits_in_list += 1
                        else:
                            decoy_hits_in_list +=1
                    fdr = (decoy_hits_in_list * self.PIT) / library_hits_in_list
                    self.fdrlist.append(fdr)
            #            print("Found fdr ", fdr, "for score ", score.score)

        #   print(fdr_at_index)
            print(len(self.fdrlist), "FDRs created")
        
    # creates a list of q values using the fdr list and sorted scores   
    
    def create_q_values(self):
        if self.qvalues:
            print("Qvalues are already created for this set")
            return False
        
        else:
            scores = [x for x in self.sorted_scores if x.type == "library"]
            for num, score in self.enumerate_reversed(scores):        
        #       print(num)
                if(num == 0):
                    qvalue = self.fdrlist[0]
        #                print("For score", score.score, " qvalue is: ", qvalue )
                    self.qvalues.append(qvalue)
                else: 
                    list_of_fdr_to_consider = self.fdrlist[ 0 : num]
                    qvalue = min(list_of_fdr_to_consider)   
        #                print("For score", score.score, " qvalue is: ", qvalue )
                    self.qvalues.append(qvalue)

            print(len(self.qvalues), " qvalues created")
            
    def clear_analysis(self):
        self.fdrlist = []
        self.qvalues = []
        

        
        