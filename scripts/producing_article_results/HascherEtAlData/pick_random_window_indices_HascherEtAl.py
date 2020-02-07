#!/usr/bin/env python3

import numpy
import numpy.random



for i in range(0,6):
    FOLDER="/path/to/NSCLC/preanalysis/chr%s"%(i)
    N_C_data=numpy.loadtxt("%s/number_of_cytosines_in_window_chr%s.txt"%(FOLDER,i))
    
    windows_with_N_C_over10=numpy.nonzero(N_C_data>=10)[0]
    
    if len(windows_with_N_C_over10)>500:
    
        random_inds=numpy.random.choice(len(windows_with_N_C_over10),500,replace=False)

        INDS_TO_BE_SAVED_ONE_BASED=windows_with_N_C_over10[numpy.sort(random_inds)]+1

        numpy.savetxt("%s/random_window_inds_N_cytosines_10_or_more_one_based_chr%s.txt"%(FOLDER,i),INDS_TO_BE_SAVED_ONE_BASED ,fmt='%d')

    else:
        print("The number of windows with N_C>=10 was smaller than 500 for chromosome %s."%(i))
    
