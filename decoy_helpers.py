from matchms.Spikes import Spikes
from custom_fragment import FragmentPeak
import random as rand

def ispeakhere(s, precursor): 
    tolerance = (precursor / 1000000) * 5
    return np.any((s.peaks.mz >= precursor - 5) & (s.peaks.mz <= precursor + 5))                                         


def masswithin5ppm(m,range):
    ishere = 0
    tolerance = (m / 1000000) * 5
    for r in range:
        if m >= r - tolerance and m <= r + tolerance:
            ishere = 1
            break
    return ishere


from custom_fragment import index, find_lt, find_le, find_gt, find_ge, find_frags_in_range


def random_sample_5_peaks(spectrum_id, all_fragments_list):
    frag_list =[x for x in all_fragments_list if x.spectrum_id == spectrum_id]
    frag_list = rand.sample(frag_list, 5)
    return frag_list


def get_spectrums_with_peak(mz, all_fragments_list_sorted):
    tolerance = (mz / 1000000) * 5
    x = FragmentPeak(mz - tolerance, 0, 'counter')
    y = FragmentPeak(mz + tolerance, 0, 'counter2')
    x = find_ge(all_fragments_list_sorted, x)
    y = find_le(all_fragments_list_sorted, y)
    
   # print("Mass for search is:", mz)
    
    first_frags = find_frags_in_range(all_fragments_list_sorted, x, y)
    
    # print("In first frags there are: ", len(first_frags), "there masses are ")   

    # for f in first_frags:
     #   print(f.mz)
        
  #  print("Now retrieving spectrums:")
        
    first_spectrum_ids = [f.spectrum_id for f in first_frags]
    
   # print("unsorted spectrum ids: ", len(first_spectrum_ids))
    
    final_spectrum_ids = list(dict.fromkeys(first_spectrum_ids))
    
    if(len(final_spectrum_ids) > 25):
        final_spectrum_ids = rand.sample(final_spectrum_ids, 25)
       # print("Too large: random sampled 25")
    
    
    #print("de-duplicated spectrum ids", len(final_spectrum_ids))
    
          
    return final_spectrum_ids



def return_random_pick(candidatefrags,decoypeaks, parentmass):
    
   # print("Number of candidates on arrival: ", len(candidatefrags)) 
    candidatefrags2 = [x for x in candidatefrags if x.mz < parentmass] 
    
    if(candidatefrags2 == None):
        candidatefrags2 = candidatefrags
    else:
        candidatefrags3 = [ x for x in candidatefrags2 if masswithin5ppm(x.mz, decoypeaks) == 0]
  
    try:
        candidateion = rand.choice(candidatefrags3)
        
    except:
        candidateion = rand.choice(candidatefrags2)

    
           
    return candidateion