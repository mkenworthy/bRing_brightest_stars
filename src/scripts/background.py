import numpy as np
import ephem
import astropy.stats as aps
import cameras

## Input: daily data
## Removes 3 sigma outliers of clouds
## Output: data without clouds
def dailycloud(datad):    
    meancloud, mediancloud, stdcloud = aps.sigma_clipped_stats(datad['clouds0'])
    datagoodloc = (datad['clouds0'] < (3*stdcloud+meancloud)) * (datad['clouds0'] > (meancloud-3*stdcloud))
    datagood = datad[datagoodloc]
    return datagood
 
## Input: daily data
## Calibrates data without 5 sigma magnitude errors
## Output: daily data without 5 sigma magnitude errors
def dailymag(data):    
    meanmag, medianmag, stdmag = aps.sigma_clipped_stats(data['emag0'])
    datagoodloc = (data['emag0'] < (5*stdmag+meanmag)) * (data['emag0'] > (meanmag-5*stdmag))
    datagood = data[datagoodloc]
    return datagood

#Input: daily data
#Removing 5 sigma outliers of sky background
#Output: daily data without 5 sigma outliers of sky background
def sky(data):
    meansky, mediansky, stdsky = aps.sigma_clipped_stats(data['sky'])    
    datagoodloc = (data['sky'] < (5*stdsky + meansky)) * (data['sky'] > (meansky-5*stdsky))
    datagood = data[datagoodloc]
    return datagood  
    
def add_sunmoon(data, camera):
    sun = ephem.Sun()
    moon = ephem.Moon()
    sunalt = data['jd'].copy()
    moonalt = data['jd'].copy()
    moonph = data['jd'].copy()
    moonmag = data['jd'].copy()
    cam = np.where(camera == np.array(cameras.cams))[0][0]
    obs = cameras.locsPE[cam]
    for iii in range(len(data)):
        obs.date = data['jd'][iii]-2415020
        moon.compute(obs)
        sun.compute(obs)
        moonalt[iii] = moon.alt*180/np.pi
        moonph[iii]  = moon.phase
        moonmag[iii] = moon.mag
        sunalt[iii]  = sun.alt*180/np.pi
    data = np.lib.recfunctions.append_fields(data, 'moonalt', moonalt)
    data = np.lib.recfunctions.append_fields(data, 'sunalt', sunalt)
    data = np.lib.recfunctions.append_fields(data, 'moonph', moonph)
    data = np.lib.recfunctions.append_fields(data, 'moonmag', moonmag)
    return data


