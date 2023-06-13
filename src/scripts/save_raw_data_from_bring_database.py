import sys
import numpy as np

# add the directories for the bring cameras
sys.path.append('/disks/web1/bring/products/')
sys.path.append('/data2/stuik/code')
import getdataALL
import getdata
import astropy.io.fits as fits
import calibration
import background

# get the star catalogue
starcat    = fits.getdata('/net/belterwijde/data2/bring/Catalogues/bringcat20191231.fits')

# Save uncalibrated data using ASCC-numbers
def save_uncalibrated(star):
        starcatloc = np.where(starcat['ascc'] == str(star))[0][0]
        starinfo = starcat[starcatloc]
        HD = starinfo['HD']
        AUW = '/data2/hoogenboom/RawData/'+str(HD)+'_AUW.fits'
        AUE = '/data2/hoogenboom/RawData/'+str(HD)+'_AUE.fits'
        SAW = '/data2/hoogenboom/RawData/'+str(HD)+'_SAW.fits'
        SAE = '/data2/hoogenboom/RawData/'+str(HD)+'_SAE.fits'

        data_AUW_daily = getdata.fromdaily(str(star), 'AUW')
        data_AUW_quarterly = getdataALL.fromquarterly(str(star), 'AUW')
        datad_AUW = background.add_sunmoon(data_AUW_daily, 'AUW')
        caldaily_AUW = calibration.get_calibration(str(star), datad_AUW, data_AUW_quarterly, 'AUW')
        fits.writeto(AUW, caldaily_AUW, overwrite = True)

        data_AUE_daily = getdata.fromdaily(str(star), 'AUE')
        data_AUE_quarterly = getdataALL.fromquarterly(str(star), 'AUE')
        datad_AUE = background.add_sunmoon(data_AUE_daily, 'AUE')
        caldaily_AUE = calibration.get_calibration(str(star), datad_AUE, data_AUE_quarterly, 'AUE')
        fits.writeto(AUE, caldaily_AUE, overwrite = True)

        data_SAW_daily = getdata.fromdaily(str(star), 'SAW')
        data_SAW_quarterly = getdataALL.fromquarterly(str(star), 'SAW')
        datad_SAW = background.add_sunmoon(data_SAW_daily, 'SAW')
        caldaily_SAW = calibration.get_calibration(str(star), datad_SAW, data_SAW_quarterly, 'SAW')
        fits.writeto(SAW, caldaily_SAW, overwrite = True)

        data_SAE_daily = getdata.fromdaily(str(star), 'SAE')
        data_SAE_quarterly = getdataALL.fromquarterly(str(star), 'SAE')
        datad_SAE = background.add_sunmoon(data_SAE_daily, 'SAE')
        caldaily_SAE = calibration.get_calibration(str(star), datad_SAE, data_SAE_quarterly, 'SAE')
        fits.writeto(SAE, caldaily_SAE, overwrite = True)

variableASCC =  [1875903, 2120290, 2198437, 2203013, 2298409, 2311128]
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
for i in nonvariableASCC:
        save_uncalibrated(i)
for i in variableASCC:
        save_uncalibrated(i)
