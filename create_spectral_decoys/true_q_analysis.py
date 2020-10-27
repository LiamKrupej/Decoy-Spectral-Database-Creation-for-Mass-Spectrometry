import sys, os
from matplotlib import pyplot as plt
from create_spectral_decoys.custom_filtering import find_chem_string, are_peaks_similar, are_spectrums_same
from create_spectral_decoys.analysis import Analysis
import numpy

# class implements Analysis class


class TrueAnalysis(Analysis):    
    def __init__(self, sorted_scores, PIT, label):
        self.sorted_scores = sorted_scores
        self.falsescores = []
        self.truescores = []
        self.fdrlist = []
        self.qvalues = []
        self.label = label
        self.PIT = PIT
        
        # elements include a list of sorted scores, total false positives and total true hits
        # FDR list, qvalues, and a label. PIT set to 0 as not used in method implementation.
        
        for s in sorted_scores:
            if(are_spectrums_same(s.query, s.library) == True):
                self.truescores.append(s.score)
            else:
                self.falsescores.append(s.score)
                
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
            for score in self.sorted_scores:
                list_to_consider = [x for x in self.sorted_scores if x.score >= score.score]
                false_hits_in_list = 0
                true_hits_in_list = 0
                for s in list_to_consider:
                    spec1 = s.query
                    spec2 = s.library
                    if(are_spectrums_same(spec1, spec2) == True):
                        true_hits_in_list += 1
                    else:
                        false_hits_in_list +=1   
                fdr = false_hits_in_list / len(list_to_consider)
                self.fdrlist.append(fdr)
        #        print("Found fdr ", fdr, "for score ", score.score)
        
        print(len(self.fdrlist), " fdrs created")
    
    # creates a list of q values using the fdr list and sorted scores
    
    def create_q_values(self):
        if self.qvalues:
            print("Qvalues are already created for this set")
            return False
        
        else:
            for num, score in self.enumerate_reversed(self.sorted_scores):    
    #            print(num)
                if(num == 0):
                    qvalue = self.fdrlist[0]
                    self.qvalues.append(qvalue) 
    #                print("For score", score.score, " qvalue is: ", qvalue )

                else: 
                    list_of_fdr_to_consider = self.fdrlist[ 0 : num]
                    qvalue = min(list_of_fdr_to_consider)
                    self.qvalues.append(qvalue) 
    #                print("For score", score.score, " qvalue is: ", qvalue)

            

        print(len(self.qvalues), " qvalues created")
        
    def clear_analysis(self):
        self.fdrlist = []
        self.qvalues = []