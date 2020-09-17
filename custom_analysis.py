from custom_filtering import find_chem_string, are_peaks_similar, are_spectrums_same
from matplotlib import pyplot as plt
import numpy


def look_for_inchi_and_precursor(query, library):
    queries_with_inchi_match = []
    queries_with_precursor_match = []
    queries_with_inchi_match_and_precursor = []

    for spec in query:
        for s in library:
            if (find_chem_string(spec) == find_chem_string(s)):
                queries_with_inchi_match.append(spec)
            if (are_peaks_similar(spec.metadata['precursor_mz'], s.metadata['precursor_mz']) == True):
                queries_with_precursor_match.append(spec)
            if (find_chem_string(spec) == find_chem_string(s)) and (
            (are_peaks_similar(spec.metadata['precursor_mz'], s.metadata['precursor_mz']) == True)):
                queries_with_inchi_match_and_precursor.append(spec)

    print(len(queries_with_inchi_match))
    print(len(queries_with_precursor_match))
    print(len(queries_with_inchi_match_and_precursor))
    type_of_match = ["Inchi_Match", "Precurs_Match", "Inchi_&_Precurs_Match"]
    num_of_matches = [len(queries_with_inchi_match), len(queries_with_precursor_match),
                      len(queries_with_inchi_match_and_precursor)]

    fig = plt.figure(figsize=(12, 7))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.bar(type_of_match, num_of_matches)
    plt.xlabel("number of peaks in spectrum")
    plt.ylabel("number of spectra in respective bin")
    plt.show()


from matchms.similarity import CosineGreedy

class CosineHit:   
    score = 0
    type = None
    query = None
    library = None
    qvalue = None
    
    def __init__(self, score, type, query, library):
        self.score = score
        self.type = type
        self.query = query
        self.library =library
        
    def __lt__(self, other):
        return self.score < other.score

    def _eq_(self, other):
        return self.score == other.score
    
    def getQuery(self):
        return self.query
    
    def getLibrary(self):
        return self.library
    
    def getQvalue(self):
        return self.qvalue
    
    def setValue(self, q):
        self.qvalue = q
        
    def setType(Self, type):
        self.type = type

def plot_hits(cosine_hits):
    
    scores = [score in cosine_hit in cosine_hits] 
    plt.figure(figsize=(12,7))
    hist = plt.hist(scores, np.arange(0,1, 0.05))
    plt.xlabel("cosine score of matches")
    plt.ylabel("number of queries in respective bin")



def return_list_cosine_scores(query, library, librarytype):
         
        if(librarytype != "library" and librarytype != "decoy"):
            print("library type paramater must be either target or decoy")
            return False
        else: 
            cosine_greedy = CosineGreedy(tolerance=0.2)
            counter = 1
            scores = []
            average_matches = 0
            milestone = 1
            
            if(librarytype == "decoy"):        
                for spec in query:
                    prelim_scores = []
                    for d in library:
                        score, n_matches = cosine_greedy(d, spec)
                        average_matches = average_matches + n_matches
                        newscore = CosineHit(score, librarytype, spec, d)
                        prelim_scores.append(newscore)          

                    prelim_scores = sorted(prelim_scores)
                    scores.append(prelim_scores[-1])
                    print(milestone)
                    milestone += 1
                    
            if(librarytype == "library"):
                for spec in query:
                    prelim_scores = []
                    for d in library:
                        if(are_peaks_similar(spec.metadata['precursor_mz'], d.metadata['precursor_mz']) == True):
                            score, n_matches = cosine_greedy(d, spec)
                            average_matches = average_matches + n_matches
                            newscore = CosineHit(score, librarytype, spec, d)
                            prelim_scores.append(newscore)   

                    prelim_scores = sorted(prelim_scores)
                    print("Scores are ")
                    for s in prelim_scores:
                        print(s.score)
                    for s in prelim_scores: 
                        if(are_spectrums_same(s.query, s.library) == True):
                            print("true score is",  s.score)
                    scores.append(prelim_scores[-1])
                    print("Score taken: ",prelim_scores[-1].score )
                    print(milestone)
                    milestone += 1

            return scores

def create_fdr_list_decoys(sorted_scores): 
    fdr_at_index = []
    for score in sorted_scores:
        list_to_consider = [x for x in sorted_scores if x.score >= score.score]
        decoy_hits_in_list = 0
        library_hits_in_list = 0
        for s in list_to_consider:     
            if(s.type == "library"):
                library_hits_in_list += 1
            else:
                decoy_hits_in_list +=1
        fdr = decoy_hits_in_list / library_hits_in_list
        if(s.type == "library"):
            fdr_at_index.append(fdr)
            print("Found fdr ", fdr, "for score ", score.score)
    
    print(fdr_at_index)
    return fdr_at_index


def create_fdr_list_for_trues(sorted_scores):   
    fdr_at_index = []
    for score in sorted_scores:
        list_to_consider = [x for x in sorted_scores if x.score >= score.score]
        false_hits_in_list = 0
        true_hits_in_list = 0
        for s in list_to_consider:
            spec1 = s.query
            spec2 = s.library
            if(are_spectrums_same(spec1, spec2) == True):
                true_hits_in_list += 1
            else:
                false_hits_in_list +=1             
        fdr = false_hits_in_list / (true_hits_in_list + false_hits_in_list)
        fdr_at_index.append(fdr)
        print("Found fdr ", fdr, "for score ", score.score)
        
    return fdr_at_index


def enumerate_reversed(seq):
    n = len(seq)
    for obj in reversed(seq):
        n -= 1
        yield n, obj
        

def create_q_values_from_fdr_list(fdrlist, scores):
    qvalues = []
    for num, score in enumerate_reversed(scores):
        if(score.type == "library"):
            print(num)
            if(num == 0):
                qvalue = fdrlist[0]
                print("For score", score.score, " qvalue is: ", qvalue )
                qvalues.append(qvalue)
            else: 
                list_of_fdr_to_consider = fdrlist[ 0 : num]
                qvalue = min(list_of_fdr_to_consider)   
                print("For score", score.score, " qvalue is: ", qvalue )
                qvalues.append(qvalue)
                     
    return qvalues


def create_q_values_from_fdr_list_true(fdrlist, scores):
    qvalues = []
    for num, score in enumerate_reversed(scores):    
        if(are_spectrums_same(score.query, score.library) == True):
            print(num)
            if(num == 0):
                qvalue = fdrlist[0]
                print("For score", score.score, " qvalue is: ", qvalue )

            else: 
                list_of_fdr_to_consider = fdrlist[ 0 : num]
                qvalue = min(list_of_fdr_to_consider)   
                print("For score", score.score, " qvalue is: ", qvalue)
                
        qvalues.append(qvalue) 
                      
    return qvalues
