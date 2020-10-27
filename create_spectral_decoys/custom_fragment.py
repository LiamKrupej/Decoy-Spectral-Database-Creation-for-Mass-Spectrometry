# a class for storing fragments with their masses and intensities


class FragmentPeak:
    def __init__(self, mz, intensity, spectrum_id):
        self.mz = mz
        self.intensity = intensity
        self.spectrum_id = spectrum_id

    def __lt__(self, other):
        return self.mz < other.mz

    def _eq_(self, other):
        return self.mz == other.mz

    def print_fragment(self):
        print("Peak with mz: ", self.mz, ", intensity: ", self.intensity, ", from spectrum: ", self.spectrum_id)

import operator
from bisect import bisect_left, bisect_right
from numpy.random import seed
from numpy.random import rand

# a set first of some standard methods from the python bisect library

def index(a, x):
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    raise ValueError

def find_lt(a, x):
    i = bisect_left(a, x)
    if i:
        return a[i-1]
    raise ValueError

def find_le(a, x):
    i = bisect_right(a, x)
    if i:
        return a[i-1]
    raise ValueError

def find_gt(a, x):
    i = bisect_right(a, x)
    if i != len(a):
        return a[i]
    raise ValueError

def find_ge(a, x):
    i = bisect_left(a, x)
    if i != len(a):
        return a[i]
    raise ValueError
    
# custom functions for bisect 
    
def find_frags_in_range(a, x, y):  
    i = bisect_left(a, x)
    j = bisect_right(a,y)
    return a[i:j]

def find_randomfrags_in_spec(fraglist, id): 
    indices = [i for i, x in enumerate(fraglist) if x.spectrum_id == id]
    frags = [fraglist[i] for i in indices]
    randomfrags = random.sample(frags, 5)
    return randomfrags