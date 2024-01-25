import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np
import paths

starcat    = fits.getdata(paths.data / 'bringcat20191231.fits.gz')

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

def openfits(filename):
  image_data = fits.getdata(filename)
  return image_data

def plotLightCurve(star):
    #load in final data
    starcatloc = np.where(starcat['ascc'] == str(star))[0][0]
    starinfo = starcat[starcatloc]
    HD = starinfo['HD']    
#    final_data = openfits('FinalData/'+str(HD)+'.fits')
    fin = str(HD)+'.fits'
    final_data = openfits(paths.data / fin)
    AUWcam = np.where(final_data['camera'] == 1)
    AUW = final_data[AUWcam]
    AUEcam = np.where(final_data['camera'] == 0)
    AUE = final_data[AUEcam]
    SAWcam = np.where(final_data['camera'] == 8)
    SAW = final_data[SAWcam] 
    SAEcam = np.where(final_data['camera'] == 7)
    SAE = final_data[SAEcam]
    
    plt.figure()
    plt.rcParams.update({'font.size': 14})

    stu = dict(alpha=0.1, rasterized=True,marker='o',mec='none',markersize=1,linestyle="")

    mj=2400000.5

    fig, axs = plt.subplots(5, figsize = (16.0, 10.0), sharex = True)
    fig.suptitle('ASCC '+ str(star)+ ', HD ' + str(HD))
    axs[0].plot(SAE['jd']-mj, SAE['mag0'], color = 'red', **stu)
    axs[0].plot(AUE['jd']-mj, AUE['mag0'], color = 'gold', **stu)
    axs[0].plot(SAW['jd']-mj, SAW['mag0'], color = 'turquoise', **stu)
    axs[0].plot(AUW['jd']-mj, AUW['mag0'], color = 'royalblue', **stu)

#    leg = axs[0].legend(['SAE', 'AUE', 'SAW', 'AUW'], bbox_to_anchor = (1, 1))

    axs[1].plot(SAE['jd']-mj, SAE['mag0'], color = 'red', **stu)
    axs[2].plot(AUE['jd']-mj, AUE['mag0'], color = 'gold', **stu)
    axs[3].plot(SAW['jd']-mj, SAW['mag0'], color = 'turquoise', **stu)
    axs[4].plot(AUW['jd']-mj, AUW['mag0'], color = 'royalblue', **stu)

    # axs[0].set(xlabel = " ", ylabel = "Magnitude")
    # axs[1].set(xlabel = " ", ylabel = "Magnitude")
    # axs[2].set(xlabel = " ", ylabel = "Magnitude")
    # axs[3].set(xlabel = " ", ylabel = "Magnitude")

    # axs[4].set(xlabel = "Time (Julian Dates)", ylabel = "Magnitude")

    fig.supylabel('Magnitude',fontsize=20)

    axs[4].set_xlabel("Time (Modified Julian Date)",fontsize=20)


    axs[0].set_title("All cameras", fontsize = 14)
    axs[1].set_title("SAE", fontsize = 14)
    axs[2].set_title("AUE", fontsize = 14)
    axs[3].set_title("SAW", fontsize = 14)
    axs[4].set_title("AUW", fontsize = 14)

    axs[0].invert_yaxis()
    axs[1].invert_yaxis()
    axs[2].invert_yaxis()
    axs[3].invert_yaxis()
    axs[4].invert_yaxis()

    fout = str(HD)+'.pdf'
    plt.tight_layout()
    plt.savefig(paths.figures / fout)
    plt.show()
    plt.close()

plotLightCurve(2311128)
