#!/usr/bin/env python
#
# Create a dataset by concatenating a series of functional images.

import sys
from bender_study import gvt

maskfile = sys.argv[1]
outfile = sys.argv[2]

# load and concatenate all input files
ds = gvt.loadcat(sys.argv[3:], maskfile)

# save the full dataset
ds.save(outfile)
print(f'Saved merged data to {outfile}.')
