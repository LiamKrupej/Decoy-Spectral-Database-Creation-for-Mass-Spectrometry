import numpy as np
import random
from matchms import Spectrum, Spikes
import sys, os
sys.path.append(r'C:\Users\User\Documents\Projects\Decoy_Spectra\create_spectral_decoys')
from custom_fragment import FragmentPeak
from custom_filtering import get_parent_peak
import time
from custom_fragment import index, find_lt, find_le, find_gt, find_ge, find_frags_in_range
import random as rand

# checks if a peak is in a spectrum by a tolerance of 5 ppm

def ispeakhere(s, precursor): 
    tolerance = (precursor / 1000000) * 5
    return np.any((s.peaks.mz >= precursor - 5) & (s.peaks.mz <= precursor + 5))                                         

# a similar function for checking any mass within any range

def masswithin5ppm(m,range):
    ishere = 0
    tolerance = (m / 1000000) * 5
    for r in range:
        if m >= r - tolerance and m <= r + tolerance:
            ishere = 1
            break
    return ishere

# random samples 5 peaks from a given spectrum 

def random_sample_5_peaks(spectrum_id, all_fragments_list):
    frag_list =[x for x in all_fragments_list if x.spectrum_id == spectrum_id]
    frag_list = rand.sample(frag_list, 5)
    return frag_list

# returns all spectrums on the list of library spectrums with a given mass peak. Since this is then used to
# to random sample five fragments from each spectrum to add to candidate fragments, output list is limited to 25.
# this is to prevent a time complexity bottleneck in the algorithm. 

def get_spectrums_with_peak(mz, all_fragments_list_sorted):
    tolerance = (mz / 1000000) * 5
    x = FragmentPeak(mz - tolerance, 0, 'counter')
    y = FragmentPeak(mz + tolerance, 0, 'counter2')
    x = find_ge(all_fragments_list_sorted, x)
    y = find_le(all_fragments_list_sorted, y)
    
    first_frags = find_frags_in_range(all_fragments_list_sorted, x, y)   
    first_spectrum_ids = [f.spectrum_id for f in first_frags]    
    final_spectrum_ids = list(dict.fromkeys(first_spectrum_ids))
    
    if(len(final_spectrum_ids) > 25):
        final_spectrum_ids = rand.sample(final_spectrum_ids, 25)    
          
    return final_spectrum_ids

# randomly selects a peak from a list of candidate fragments. Ideally, peak mz is less than parent mass and not within 5ppm
# of another peak. Where no candidate frags are less than parentmass, this criterion is dropped. Where all candidate frags are
# within 5ppm of another, this criterion is dropped. 

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

# a class for creating both naive and spectrum-based decoys

class DecoyFactory:   
    
    # object is passed a list of target library spectrums in constructor
    
    def __init__(self, spectrums):
        self.target_spectrums = spectrums
        self.all_fragments_list = []
        self.n_decoys = []
        self.sb_decoys = []
        self.sb_rejects = []
        
        
     # a method for stripping away peaks of all the target library spectra, turning them to FragmentPeak objects
    # and storing them as part of a sorted fragment list
        
    
    def strip_fragments(self):
        
        if self.all_fragments_list:
            print("spectra already fragmented")
            return False
        print("Beginning fragmentation process....")  
        for spec in self.target_spectrums:
            i = 0;
            for p in spec.peaks.mz:
                frag = FragmentPeak(spec.peaks.mz[i], spec.peaks.intensities[i], spec.get('spectrum_id'))
                i += 1
                self.all_fragments_list.append(frag)                    
        self.all_fragments_list = sorted(self.all_fragments_list)
        print("Number of fragments made ", len( self.all_fragments_list))
        
    # creates a naive decoy using target library fragment list and passed spectrum
    
    def create_naive_decoy(self, s):
        # print(get_parent_peak(s))
        decoy_mz = np.array([get_parent_peak(s)[0]])
        decoy_intensity = np.array([get_parent_peak(s)[1]])  
        metadata = {'precursor_mz': decoy_mz[0]} 
        peaks_in_target = len(s.peaks.mz)
        random_spectrums = random.sample(self.target_spectrums, peaks_in_target - 1)
        for spec in random_spectrums:
            randommass =  random.choice(spec.peaks.mz)
            index = np.where(spec.peaks.mz == randommass)
            randomintensity = spec.peaks.intensities[index]
            decoy_mz = np.append(decoy_mz, [randommass])
            decoy_intensity = np.append(decoy_intensity, [randomintensity])
        decoy_mz = np.asarray(decoy_mz, dtype=float) 
        decoy_intensity = np.asarray(decoy_intensity, dtype=float) 
        inds  = decoy_mz.argsort()
        sorted_intensities = decoy_intensity[inds]
        sorted_mzs = decoy_mz[inds]
        decoy = Spectrum(sorted_mzs, sorted_intensities, metadata)
        return decoy
    
    # creates a naive deoy spectrum for every spectrum in the target library
                                    
    def create_naive_decoys_list(self):                
        if self.n_decoys:
            print("Naive Decoys already created")
            return False                  
        decoysprocessed = 1
        naive_decoy_spectrums = []

        for spec in self.target_spectrums:
            try:
                s = self.create_naive_decoy(spec)
                self.n_decoys.append(s)
            except:
                print("failed to create decoy ", (decoysprocessed -1))
            decoysprocessed += 1
            milestones = [10, 100, 500, 1000, 1500, 2000]
            for m in milestones:
                    if(decoysprocessed == m):
                        print(m, " naive decoys processed")
                     
        print(len(self.n_decoys), " naive decoys created.")
        
        
        
     # creates a spectrum-based decoy
     
            
    def create_spectrum_based_decoy(self, s):
       # print("This spectrum has: ", len(s.peaks.mz), " peaks.")
        parentmass = get_parent_peak(s)[0]
        parentintensity = get_parent_peak(s)[1]
        metadata = {'precursor_mz': parentmass} 
        decoy_mz = np.array([parentmass])
        decoy_intensities = np.array([parentintensity])
    #    print("Parent peak equals: ", parentmass, "m/z, with intensity: ", parentintensity)
        peaks_in_target = len(s.peaks.mz)  
        candidate_fragments_list = []
        mass_for_loop_seeding = parentmass.copy()

        while(len(decoy_mz) < len(s.peaks.mz)):

            id_list = get_spectrums_with_peak(mass_for_loop_seeding, self.all_fragments_list)

            for id in id_list:
                random_peaks = random_sample_5_peaks(id, self.all_fragments_list)
                candidate_fragments_list.extend(random_peaks)

    #        print("Length of candidate frags list: ", len(candidate_fragments_list))
            drawn_ion = return_random_pick(candidate_fragments_list, decoy_mz, parentmass)

    #        print("Drew randomly:", drawn_ion.mz)
            decoy_mz = np.append(decoy_mz, drawn_ion.mz)
            decoy_intensities = np.append(decoy_intensities, drawn_ion.intensity)


    #        print("Added peak with mass ", drawn_ion.mz, "and intensity ", drawn_ion.intensity)
    #        print("Decoy mz is length ", len(decoy_mz))

            mass_for_loop_seeding = drawn_ion.mz       

    #    print("Decoy masses has this number: ", len(decoy_mz))
    #    print("Decoy intensities has this number: ", len(decoy_intensities)) 

        decoy_mz = np.asarray(decoy_mz, dtype=float) 
        decoy_intensities = np.asarray(decoy_intensities, dtype=float) 
        inds  = decoy_mz.argsort()
        sorted_intensities = decoy_intensities[inds]
        sorted_masses = decoy_mz[inds]
        decoy = Spectrum(sorted_masses, sorted_intensities, metadata) 


        return decoy
    
    
    # creates one spectrum based decoy for every spectrum in the target library

    def create_spectrum_based_decoy_list(self):
        if self.sb_decoys:
            print("Spectrum-Based Decoys already created")
            return False  
        else:
            
            decoysprocessed = 1
            times_taken = []
      
            for spec in self.target_spectrums:
                try:
                    s = self.create_spectrum_based_decoy(spec)
                    self.sb_decoys.append(s)
                except:  
                    self.sb_rejects.append(spec)
                decoysprocessed += 1
                milestones = [10, 100, 500, 1000, 1500, 2000]
                for m in milestones:
                    if(decoysprocessed == m):
                        print(m, " decoys processed")
                        print(len(self.sb_decoys), " SB decoys created")
            print("Created: ", len(self.sb_decoys), ". Failed to create: ", len(self.sb_rejects))
            
    def clear_naive_decoys(self):
        self.n_decoys = []
        
    def clear_sb_decoys(self):
        self.sb_decoys = []
        self.sb_rejects = []
    
    def clear_fragments(self):
        self.fragments = []
        
        
    
                  
                  