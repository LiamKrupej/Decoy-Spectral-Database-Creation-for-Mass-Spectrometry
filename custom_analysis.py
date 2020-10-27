from create_spectral_decoys.custom_filtering import find_chem_string, are_peaks_similar, are_spectrums_same
from matplotlib import pyplot as plt
import numpy

# takes a list of query spectra and library spectra and returns numbers of inchi matches and precursor matches in the form of 
# bar chart

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

# a class to store hits between query and library spectrum lists

class CosineHit:   
    def __init__(self, score, type, query, library):
        self.score = score
        self.type = type
        self.query = query
        self.library =library
        
    def __lt__(self, other):
        return self.score < other.score

    def _eq_(self, other):
        return self.score == other.score

# takes query list and library list and scores for similarity and stores highest score per query on a CosineHit list

def return_list_cosine_scores(query, library, type):
         
    if(type != "library" and type != "decoy"):
        print("library type parameter must be either library or decoy")
        return False
    else: 
        cosine_greedy = CosineGreedy(tolerance=0.2)
        counter = 1
        scores = []
        average_matches = 0
        milestone = 1

        if(type == "decoy"):        
            for spec in query:
                prelim_scores = []
                for d in library:
                    score, n_matches = cosine_greedy(d, spec)
                    average_matches = average_matches + n_matches
                    newscore = CosineHit(score, type, spec, d)
                    prelim_scores.append(newscore)          

                prelim_scores = sorted(prelim_scores)
                scores.append(prelim_scores[-1])


        if(type == "library"):
            for spec in query:
                prelim_scores = []
                for d in library:
                    if(are_peaks_similar(spec.metadata['precursor_mz'], d.metadata['precursor_mz']) == True):
                        score, n_matches = cosine_greedy(d, spec)
                        average_matches = average_matches + n_matches
                        newscore = CosineHit(score, type, spec, d)
                        prelim_scores.append(newscore)  
                    else:
                        newscore = CosineHit(0, type, spec, d)

                prelim_scores = sorted(prelim_scores)
#                    print("Scores are ")
#                    for s in prelim_scores:
#                        print(s.score)
#                    for s in prelim_scores: 
#                        if(are_spectrums_same(s.query, s.library) == True):
#                            print("true score is",  s.score)
                scores.append(prelim_scores[-1])
#                    print("Score taken: ",prelim_scores[-1].score )



        return scores