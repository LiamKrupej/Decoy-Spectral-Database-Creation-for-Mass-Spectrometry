from abc import ABC, abstractmethod

# abstract class for analysis

class Analysis(ABC):
    
    # each class must have a set of scores, a PIT and a label
    
    def __init__(self, sortedscores, PIT, label):
        self.sortedscores = sortedscores
        self.PIT = PIT
        self.label = label   
        
        super().__init()
     
    # enumerate_reversed is for reversing sequences of objects, used in FDR and q value methods
    
    def enumerate_reversed(self, seq):
        pass
     
    # create fdr list from a given list of scores
    
    def create_fdr_list(self): 
        pass
    
    # create q values from a given list of scores
    
    def create_qvalues(self):
        pass
    
    

