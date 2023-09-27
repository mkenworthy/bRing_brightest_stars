import sys
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from numpy.lib.recfunctions import stack_arrays
import calibration
import detrendsteps


sys.path.append('/disks/web1/bring/products/')
starcat    = fits.getdata('/net/belterwijde/data2/bring/Catalogues/bringcat20191231.fits')

nonvariableASCC = [1721367, 1722082, 1728196, 1728326, 1730728, 1748375, 1762319, 1776717, 1781002, 1790528, 1821588, 1822707, 1831780, 1839751, 1862167, 
                   1862828, 1865291, 1867878, 1871021, 1873533, 1875934, 1880579, 1880898, 1881897, 1882798, 1887460, 1901232, 1905964, 1905967, 1907532, 
                   1909803, 1911993, 1923547, 1925996, 1927667, 1936013, 1938697, 1943462, 1959119, 1959126, 1959843, 1959875, 1962822, 1964703, 1967127,  
                   1967246, 1968585, 1974877, 1976808, 1980004, 1981329, 1992514, 1992594, 2000576, 2003304, 2004703, 2005323, 2021081, 2023562, 2027459, 
                   2028223, 2030506, 2040152, 2041127, 2050494, 2050813, 2056697, 2058847, 2061037, 2064147, 2079242, 2086192, 2094623, 2098110, 2099216, 
                   2100006, 2104402, 2110269, 2111805, 2113099, 2117254, 2119926, 2129944, 2141153, 2146573, 2146914, 2148632, 2154691, 2162355, 2182568, 
                   2195224, 2199019, 2212858, 2217847, 2218411, 2219507, 2229331, 2234363, 2237409, 2246690, 2248482, 2250231, 2268985, 2274286, 2275069, 
                   2277091, 2277100, 2286020, 2286969, 2290853, 2293937, 2295857, 2296439, 2300928, 2307136, 2308342, 2315295, 2317331, 2318884, 2325848, 
                   2333058, 2333718, 2345078, 2348878, 2349085, 2354207, 2356113, 2361416, 2362160, 2370345, 2381455, 2383925, 2385516, 2386073, 2387890, 
                   2395368, 2399317, 2399980, 2410740, 2417375, 2424942, 2425895, 2430714, 2438259, 2445952, 2455523, 2456263, 2472748, 2481765, 2493279, 
                   2497346]

variableASCC =  [1875903, 2120290, 2198437, 2203013, 2298409, 2311128]
periodVariableASCC = [1.44626907, 45.1503, 1.6697671, 2.943, 9.8426, 35.53584]

## After calibrating and detrending data, plot the lightcurves
def saveLightCurve(star, variability = False, period = 0):
    starcatloc = np.where(starcat['ascc'] == str(star))[0][0]
    starinfo = starcat[starcatloc]
    HD = starinfo['HD']
    
    dataAUW = calibration.calibrate(star, 'AUW')
    dataAUE = calibration.calibrate(star, 'AUE')
    dataSAW = calibration.calibrate(star, 'SAW')
    dataSAE = calibration.calibrate(star, 'SAE')
    
    #Detrending data
    if variability == False:
        #Detrending not variable stars
        detrendedStars = detrend(dataAUW, dataAUE, dataSAW, dataSAE)
        
    else:
        #Detrending variable stars
        detrendedStars = detrendVariable(dataAUW, dataAUE, dataSAW, dataSAE, period)

    #save final files
    finaldata = detrendedStars
    newarray = np.recarray(len(finaldata), dtype = finaldata.dtype)
    columns = newarray.dtype
    for n in columns.names:
       newarray[n] = finaldata[n]
    #  save newarray using fits.writeto
    fits.writeto('/data2/hoogenboom/FinalData/'+str(HD)+'.fits', newarray, overwrite = True)
    return detrendedStars

## Detrending the non-variable stars
## First we calculate the offset and next the offset is used to detrend the data
def detrend(dataAUW, dataAUE, dataSAW, dataSAE):
    data = [dataAUW, dataAUE, dataSAW, dataSAE] 
    detrendData = []

    for i in data:
        offsetLST = detrendsteps.lstday(i, np.zeros(len(i)))[1]
        offsetJD = detrendsteps.jdday(i, offsetLST)[1]
        offsetLSTJD = offsetLST + offsetJD
        offsetMoon = detrendsteps.moon(i, offsetLSTJD)[1]
        offset = offsetMoon + offsetLSTJD

        offsetLST2 = detrendsteps.lstday(i, offset)[1]
        offsetLST2 += offset
        offsetJD2 = detrendsteps.jdday(i, offsetLST2)[1]
        offsetLSTJD2 = offsetLST2 + offsetJD2
        detrendFinal = detrendsteps.moon(i, offsetLSTJD2)[0]
        detrendData.append(detrendFinal)
    
    detrendstack = stack_arrays((detrendData[0], detrendData[1], detrendData[2], detrendData[3]), asrecarray=True, usemask=False)
    return detrendstack

## Detrending the variable stars
## First we calculate the offset and next the offset is used to detrend the data
def detrendVariable(dataAUW, dataAUE, dataSAW, dataSAE, period):
    data = [dataAUW, dataAUE, dataSAW, dataSAE]
    detrendData = []

    for i in data:
        offsetVar = detrendsteps.variable(i, np.zeros(len(i)), period)[1]
        offsetLST = detrendsteps.lstday(i, offsetVar)[1]
        offsetLSTVar = offsetVar + offsetLST
        offsetJD = detrendsteps.jdday(i, offsetLSTVar)[1]
        offsetLSTJD = offsetLST + offsetJD
        offsetMoon = detrendsteps.moon(i, offsetLSTJD)[1]
        offset = offsetMoon + offsetLSTJD

        offsetLST2 = detrendsteps.lstday(i, offset)[1]
        offsetLST2 += offset
        offsetJD2 = detrendsteps.jdday(i, offsetLST2)[1]
        offsetLSTJD2 = offsetLST2 + offsetJD2
        detrendFinal = detrendsteps.moon(i, offsetLSTJD2)[0]
        detrendFinal['mag0'] += offsetVar
        detrendData.append(detrendFinal)

    detrendstack = stack_arrays((detrendData[0], detrendData[1], detrendData[2], detrendData[3]), asrecarray=True, usemask=False)
    return detrendstack
    
#for i in nonvariableASCC:
#    saveLightCurve(i)
#    print(i)

#for i in range(len(variableASCC)):
#    saveLightCurve(variableASCC[i], True, periodVariableASCC[i])

data = saveLightCurve(variableASCC[5], True, periodVariableASCC[5])


