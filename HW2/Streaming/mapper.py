#!/usr/bin/env python


import sys
import math

# input from STDIN
for line in sys.stdin:
    # remove leading and trailing whitespace
    line = line.strip()
    
    # split the line into x and y 
    xy = line.split()
    
    # transform x and y to floats:
    x = float(xy[0])
    y = float(xy[1])
    
    # find the lower bin boundaries for x and y; then transform to strings:
    x_lo = str(math.floor(x*10)/10)
    y_lo = str(math.floor(y*10)/10)

    # create the key, using comma separation. Since there is a 1-to-1 relationsgip between x_lo and x_hi (and similarly between y_lo and y_hi), there's no need to make the key excessively long by including both.  We'll infer x_hi from x_lo and y_hi from y_lo in the reduce step 
    xy_lo = x_lo+','+y_lo
    
    # write the tab separated (key, value) pair to STDOUT.  Just use 1 as the count, since we haven't yet aggregated the data in any way.
    print '%s\t%s' % (xy_lo,1)

