import numpy


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


def find_chem_string(spec):
    s = spec.metadata['inchi']
    if (s == "") or (s is None):
        return None
    else:
        string_in_slashes = s.split('/')[1].strip()
        return string_in_slashes


def are_peaks_similar(mass1, mass2):
    tolerance = 1
    if (-1 <= (mass1 - mass2) <= 1):
        return True
    else:
        return False


def are_spectrums_same(spec1, spec2):
    chemstring1 = find_chem_string(spec1)
    chemstring2 = find_chem_string(spec2)
    peak1 = spec1.metadata['precursor_mz']
    peak2 = spec2.metadata['precursor_mz']

    if ((chemstring1 == chemstring2) and (are_peaks_similar(peak1, peak2)) == True):
        return True
    else:
        return False


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
