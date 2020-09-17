from collections import defaultdict
from custom_fragment import FragmentPeak
from custom_fragment import find_frags_in_range
from matchms import Spikes, Spectrum
import numpy as np
from matchms.filtering import normalize_intensities 
from custom_filtering import are_peaks_similar
from adduct_calculator.adduct_rules import AdductTransformer
from molmass import Formula

t  = AdductTransformer()


class AdductAndMassLibrary:
    chemical_formulas = []
    adducts = []
    t = AdductTransformer()
    masslibraries = {}
    
    def __init__(self, adducts, chemical_formulas):
        self.adducts = adducts
        self.chemical_formulas = chemical_formulas
        masslibraries = {}
        for a in adducts: 
            masslibraries[a] = []
        for a in adducts:
            print("...", a)
            for mol in chemical_formulas:
                f = Formula(mol)
                mass = t.mass2ion(f.mass, a)
                frag = FragmentPeak(mass, None, None)
                masslibraries[a].append(frag)
        for a in adducts:
            masslibraries[a] = sorted(masslibraries[a])
        self.masslibraries = masslibraries
    
    def is_mass_here(self, mass, adduct):
        ishere = False
        lib = self.masslibraries[adduct]  
        tolerance = 0.002
        x = FragmentPeak(mass-tolerance, None, None)
        y = FragmentPeak(mass + tolerance, None, None)
        similarpeaks = find_frags_in_range(lib, x, y)
       # if(len(similarpeaks) > 0):
           # print("Similar peaks:")
           # for p in similarpeaks:
           #     print(p.mz)
           # print("Similar peaks ended:")
        if similarpeaks:
            ishere = True
        return ishere
    
    def turn_to_fragments(self, spectrum):
        frags = []
        for index, mass in enumerate(spectrum.peaks.mz):
            intensity = spectrum.peaks.intensities[index]
            f = FragmentPeak(mass, intensity, None)
            frags.append(f)  
        return frags
    
    def filter_out_noisy_peaks(self, spectrum):
        adt = spectrum.metadata['adduct']
        frags = self.turn_to_fragments(spectrum)
        mzlist = []
        intensitylist = []
        for frag in frags:
            if(self.is_mass_here(frag.mz, adt) != False) or (are_peaks_similar(frag.mz, spectrum.metadata['precursor_mz'])):
                mzlist.append(frag.mz)
                intensitylist.append(frag.intensity)
        return mzlist, intensitylist
        
     
    def return_noise_filtered_spectrum(self, spectrum):
       # spectrum.plot()
        md = spectrum.metadata.copy()
        print("Number of peaks on entry:", len(spectrum.peaks))
        mzlist, intensitylist = self.filter_out_noisy_peaks(spectrum)
        mzlist = np.asarray(mzlist, dtype=float)
        intensitylist = np.asarray(intensitylist, dtype=float)
        spec = Spectrum(mzlist, intensitylist, md)
        spec = normalize_intensities(spec)
        print("number on exit", len(spec.peaks.mz))
       # spec.plot()
        return spec
                
                         
        print("Number of peaks on exit:", len(spectrum.peaks)) 