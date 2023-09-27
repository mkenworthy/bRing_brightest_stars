import sys
import numpy as np
import astropy.stats as aps

sys.path.append('/disks/web1/bring/products/')
sys.path.append('/data2/stuik/code')

import getdataALL

#Input: ASCC-number of the star, camera
#Receives quarterly data for a star and calibrates the data using the sunaltitude, emag0 and clouds
#Output: quarterly data without clouds and good range of sunaltitude and emag0
def quarterly(star, camera = 'AUW'):
    dataq = getdataALL.fromquarterly(str(star), camera)
    # data_goodlocq = np.where((dataq['sunalt'] < -18)) #Good data locations using emag0 and sunalt
    #data_goodq = dataq # np.array(dataq[data_goodlocq]) #Putting good data of emag0 and sunalt in an array
    meanecloud, medianecloud, stdecloud = aps.sigma_clipped_stats(dataq['eclouds0'])
    data_filter = np.where((dataq['eclouds0'] < (3*stdecloud+meanecloud)) & (dataq['eclouds0'] > (meanecloud-3*stdecloud))) #Cloud error in 3 sigma range
    data = np.array(dataq[data_filter]) #Putting good data in array
    return data

