print("Loading Test Data!")


import os
import sys
from random import random
import gensim
import numpy as np
import pandas as pd
import custom_filtering
from matplotlib import pyplot as plt
from matchms import Scores, Spectrum
from matchms.importing import load_from_json
from matchms.filtering import normalize_intensities
from matchms.filtering import require_minimum_number_of_peaks
from matchms.filtering import select_by_mz
from matchms.filtering import select_by_relative_intensity
from matchms.filtering import reduce_to_number_of_peaks
from matchms.filtering import add_losses
from matchms.filtering import reduce_to_number_of_peaks
import random
from decoy_factory import random_sample_5_peaks, get_spectrums_with_peak, return_random_pick
from custom_fragment import FragmentPeak
from decoy_factory import DecoyFactory



ROOT = os.path.dirname(os.getcwd())
#path_data = os.path.join(ROOT, 'data')
path_data = 'C:\\Users\\User\\Data'
sys.path.insert(0, ROOT)

from matchms.importing import load_from_json

filename = os.path.join(path_data,'gnps_positive_ionmode_cleaned_by_matchms_and_lookups.json')
spectrums = load_from_json(filename)

print("number of spectra:", len(spectrums))


print("Testing begun!")

               


