import numpy as np
import astropy.stats as aps

synodic = 29.530588853

#Input: 
#data - calibrated data of a star
#offset - should be an empty array the first time you use this definition. The other time you use the offset given by def moon added to the offset of lstday and jdday.
#Output: detrended dataset using lstseq and an offset array
def lstday(data, offset = []):
    
    datac = data.copy()
    datau = np.unique((datac['lstseq']%13500)) #Make an array using the unique values in lstseq%13500. This means that every time data has been taken is selected one time (there are 13500 data measurements in one day). 
    offset_lst = np.zeros(len(datac)) #Array which has the same length as the dataset
    datac['mag0'] = datac['mag0'] - offset #Removing the offset array from the values of the dataset
    jd = [2457755, 2458120, 2458485, 2458850, 2459216, 2459581, 2459946, 2460311] #1 Jan 2017, 1 Jan 2018, 1 Jan 2019, 1 Jan 2020, 1 Jan 2021, 1 Jan 2022, 1 Jan 2023, 1 Jan 2024
    for i in (datau):    
        for j in range(len(jd) -1): #Calculating the mean for every year
            true = np.where(((datac['lstseq']%13500) == i) & (datac['jd'] < jd[j+1]) & (datac['jd'] >= jd[j])) #Selecting as a boolian in the data on which integers the value of i is situated and make sure the mean is calculated in one specific year
            image = (datac)[true] #Make an array where every lstseq%13500 is the same
            mean, md, RMS_raw = aps.sigma_clipped_stats(image['mag0']) #Calculate mean of every mag0 of i
            image['mag0'] = image['mag0'] - mean #Subtract the mean value of mag0 of every mag0 of i
            np.put(datac['mag0'], true, image['mag0']) #Replace the values of mag0 on i by mag0-mean 
            np.put(offset_lst, true, mean) #Replace the empty offset array by the mean of mag0
    return datac, offset_lst


#Input: 
#data - calibrated data of a star and detrended by at least lstseq
#offset - use the offset given by at least def lstday
#Output: detrended dataset using lstseq+jd and an offset array
def jdday(data, offset):
    datac = data.copy()
    day = 500 #Decide how many datapoint you would like to have in one day
    dayidx = np.floor((datac['jd']%1) * day) #Selecting the data for one day and multiply this by day so that you have the right amount of datapoint you want in one day. Round these values down.
    offset_jd = np.zeros(len(datac))
    datac['mag0'] = datac['mag0'] - offset
    for i in (range(day)):
        loc = np.where((dayidx) == i)[0]
        image = (datac)[loc]
        mean, md, RMS_raw = aps.sigma_clipped_stats(image['mag0'])
        image['mag0'] = image['mag0'] - mean
        np.put(datac['mag0'], loc, image['mag0'])
        np.put(offset_jd, loc, mean)
    return datac, offset_jd

#Input: 
#data - calibrated data of a star
#offset - use the offset given by at least def lstday + def jdday
#Output: detrended dataset using lstseq+jd+moon and an offset array
def moon(data, offset):
    datac = data.copy()
    day = 25
    dayidx = np.floor((datac['jd'] % synodic)/synodic * day) #Selecting the datapoints so you have a period of one mooncicle
    offset_moon = np.zeros(len(datac))
    datac['mag0'] = datac['mag0'] - offset
    for i in (range(day)):
        loc = np.where((dayidx) == i)[0]
        image = (datac)[loc]
        mean, md, RMS_raw = aps.sigma_clipped_stats(image['mag0'])
        image['mag0'] = image['mag0'] - mean
        np.put(datac['mag0'], loc, image['mag0'])
        np.put(offset_moon, loc, mean)
    return datac, offset_moon

def variable(data, offset, period = 9.8426, dag = 50):
    datac = data.copy()
    dagidx = np.floor((datac['jd'] % period)/period * dag)
    offset_variable = np.zeros(len(datac))
    datac['mag0'] = datac['mag0'] - offset
    for i in (range(dag)):
        loc = np.where((dagidx) == i)[0]
        image = (datac)[loc]
        mean, md, RMS_raw = aps.sigma_clipped_stats(image['mag0'])
        image['mag0'] = image['mag0'] - mean
        np.put(datac['mag0'], loc, image['mag0'])
        np.put(offset_variable, loc, mean)
    return datac, offset_variable

