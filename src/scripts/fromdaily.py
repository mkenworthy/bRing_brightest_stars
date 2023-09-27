import sys
import numpy as np

#sys.path.append('/disks/web1/bring/products/')
#sys.path.append('/data2/stuik/code')

import getdata

#Input: ASCC-number of the star, camera
#Receives the daily data for a star and calibrates the data using short exposures, pflag and aflag
#Output: daily data without flags and long exposures
def daily(star, camera = 'AUW'): 
    datad = getdata.fromdaily(str(star), camera)
    datashortloc = np.where((datad['lstseq'] %2) == 1) #Gives short exposures locations
    datashort = np.array(datad[datashortloc]) #Gives array of short exposures
    data_goodloc = np.where((datad['pflag'] == 0) & (datad['aflag'] <= 1)) #Good data locations without pflag, aflag and emag0
    data_good = np.array(datad[data_goodloc]) #Putting good data of pflag, aflag and emag0 in an array
    datashort_goodloc = np.where((datashort['pflag'] == 0) & (datashort['aflag'] <= 1)) #Good data locations without pflag, aflag and emag0 for short exposures
    data_d = np.array(datashort[datashort_goodloc]) #Putting good data of pflag, aflag and emag0 of short exposures in an array
    return data_d


