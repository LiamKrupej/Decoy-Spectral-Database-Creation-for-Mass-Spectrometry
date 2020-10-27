import numpy

# a function that takes the precursor mz from the metdata and checks if it is represented as 
# the parent peak in the spectrum's spikes (within 5 ppm)and returns it

def get_parent_peak(s):
    precursor = s.metadata['precursor_mz']
    tolerance = (precursor / 1000000) * 5
    condition = (s.peaks.mz >= precursor - tolerance) & (s.peaks.mz <= precursor + tolerance)
    mz_candidates = s.peaks.mz[condition]
    intensities_candidates = s.peaks.intensities[condition]
    if len(mz_candidates) is 0:
        return None
    else:
        max_index = numpy.argmax(intensities_candidates)
        return mz_candidates[max_index], intensities_candidates[max_index]

# return inchi from inchi key, none if not present

def find_chem_string(spec):
    s = spec.metadata['inchi']
    if (s == "") or (s is None):
        return None
    else:
        string_in_slashes = s.split('/')[1].strip()
        return string_in_slashes

# checks if peaks are similar within a larger mass tolerance of 1 Da

def are_peaks_similar(mass1, mass2):
    tolerance = 1
    if (-1 <= (mass1 - mass2) <= 1):
        return True
    else:
        return False

# checks if two spectrums represent the same molecule by matching both inchi and precursor mass

def are_spectrums_same(spec1, spec2):
    chemstring1 = find_chem_string(spec1)
    chemstring2 = find_chem_string(spec2)
    peak1 = spec1.metadata['precursor_mz']
    peak2 = spec2.metadata['precursor_mz']

    if ((chemstring1 == chemstring2) and (are_peaks_similar(peak1, peak2)) == True):
        return True
    else:
        return False

# checks if a molecule is in a list of spectrums    

def is_molecule_here(spec, spec_list):
    ishere = 0
    mcl = find_chem_string(spec)
    for s in spec_list:
        lib_mcl = find_chem_string(s)
        mass1 = spec.metadata['precursor_mz']
        mass2 = s.metadata['precursor_mz']
        if(mcl == lib_mcl) and (are_peaks_similar(mass1, mass2) == True):
            ishere = 1
    return ishere

# takes a list of spectrums and returns only those with a parent peak / precursor match

def return_spectrums_with_parent_peaks(spectrums):
    specs = spectrums
    spectrumswithpeak = []
    for s in specs:
        if(get_parent_peak(s) != None):
            spectrumswithpeak.append(s)
    return spectrumswithpeak

# takes a list of spectrums and returns only those with an inchi

def return_spectrums_with_inchi(spectrums):
    spectrumswithinchi = []
    for spec in spectrums:
        if(find_chem_string(spec) is not None):
            spectrumswithinchi.append(spec)
    return spectrumswithinchi

# takes a list of spectrums and returns only those with an adduct

def return_spectrums_with_adduct(spectrums):
    spetrumswithadduct = []  
    for spec in spectrums:
        if(isinstance(spec.metadata['adduct'], str)) & (spec.metadata['adduct'] != ""):
            spetrumswithadduct.append(spec)
    return spetrumswithadduct

# takes a list of spectrums and a list of adducts and returns only spectrums with those adducts

def return_spectrums_with_given_adducts(spectrums, positive_adducts):
    spectrumswithposadduct = []
    for spec in spectrums:
        if(spec.metadata['adduct'] in positive_adducts):
            spectrumswithposadduct.append(spec)
    return spectrumswithposadduct


