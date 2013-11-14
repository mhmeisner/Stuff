#!/usr/bin/env python

from operator import itemgetter
import sys

# set starting values:
current_x_lo = None
current_x_hi = None
current_y_lo = None
current_y_hi = None
current_count = 0

# read lines from STDIN (the output from the map step)
for line in sys.stdin:
    # remove whitespace
    line = line.strip()
    
    #split into key and value (which were tab separated)
    xy_lo,count = line.split('\t')
    
    # split into x_lo and y_lo, which were strung together with a comma to make up the key
    x_lo = xy_lo.split(',')[0]
    y_lo = xy_lo.split(',')[1]
    
    # try converting the count to an integer...should be no problem as they all should be 1. 
    try:
        count = int(count)
    except ValueError:
		# count was not a number, so silently ignore/discard this line
        continue

    # check and see if the bin boundaries are the same as the previous line. If so, increase the count and don't print anything. Again, since there's a 1-to-1 relationship between x_lo and x_hi, we only need to check either the lower bounds of the bins (as I've chosen) or the upper bounds, for both x and y.
    if current_x_lo == x_lo and current_y_lo == y_lo:
        current_count += count
    else: # we know we've moved onto a new bin, so print out the boundaries of the previous bin and the total count of points in that bin:
        if current_x_lo and current_y_lo: # makes sure that x_lo and y_lo aren't NONE, therefore handling the first line properly, where current_x_lo and current_y_lo will not have been updated yet...we don't have a previous bin yet so we can't have moved onto a new one!  
            # convert lower bounds to floats, add .1 to get the upper bound, and return them to strings:
            current_x_hi = str(float(x_lo) + 0.1)
            current_y_hi = str(float(y_lo) + 0.1)
        # print the output in the desired comma-separated (x_lo,x_hi,y_lo,y_hi,count) format: 
            print '%s,%s,%s,%s,%s' % (current_x_lo,current_x_hi,current_y_lo,current_y_hi,current_count)
        # also update the current values of x_lo and y_lo to reflect the new boundaries of the bin in which the point printed on the current line lies in.  Reset count to whatever the count was for this line (always 1 in this case). This is outside of the previous if statement since we want to set the current values for the first time when we read in the first line 
        current_count = count
        current_x_lo = x_lo
        current_y_lo = y_lo


# output the last bin count:
# convert lower bounds to floats, add .1 to get the upper bound, and return them to strings:
current_x_hi = str(float(x_lo) + 0.1)
current_y_hi = str(float(y_lo) + 0.1)
if current_x_lo == x_lo and current_y_lo == y_lo:
    print '%s,%s,%s,%s,%s' % (current_x_lo, current_x_hi,current_y_lo,current_y_hi,current_count)

