#!/usr/bin/env python
#
# Calculate average of residual parameters.

import sys
import numpy as np

acf = np.loadtxt(sys.argv[1])
acf_mean = np.mean(acf, 0)
s = '{:.6f} {:.6f} {:.6f}'.format(acf_mean[0],acf_mean[1],acf_mean[2])
print(s)
