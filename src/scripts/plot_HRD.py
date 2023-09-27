#plot HR diagram
#uses "starcat" which should contain the same information as "gaiacat"
import paths

import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as fits
from matplotlib import colors
plt.rc('text', usetex=True)

# using histo code from https://vlas.dev/post/gaia-dr2-hrd/

starcat    = fits.getdata(paths.data / 'bringcat20191231.fits.gz')
gaiacat    = fits.getdata(paths.data / 'gaiadr2.fits.gz')

stars  = [2111805, 2348878, 2199019, 2345078, 2333718, 2250231, 1880898, 2248482, 2386073, 2098110, 2023562, 1980004, 2417375, 2286020, 2120290, 2212858, 1865291, 2099216, 1936013, 1927667, 2217847, 1781002, 2154691, 2061037, 1962822, 1881897, 1905967, 2050813, 2218411, 2056697, 2146573, 1964703, 1880579, 2399317, 2041127, 2318884, 1862167, 2246690, 2430714, 1968585, 2277091, 2354207, 2079242, 2293937, 2370345, 2410740, 2113099, 1875903, 2387890, 1901232, 1981329, 1728196, 2094623, 1822707, 2275069, 2325848, 1887460, 2219507, 2349085, 1882798, 1967127, 1909803, 2203013, 1862828, 1790528, 2300928, 1923547, 2438259, 2399980, 2455523, 2277100, 2004703, 1976808, 2296439, 2317331, 2315295, 1967246, 1959119, 2162355, 1873533, 2425895, 1907532, 2117254, 1959126, 2195224, 2086192, 1730728, 2424942, 2129944, 1762319, 2058847, 2104402, 1721367, 1875934, 1871021, 1938697, 2362160, 2119926, 2361416, 2333058, 1974877, 1831780, 2395368, 2268985, 2286969, 2182568, 1748375, 1925996, 2497346, 2311128, 2030506, 2445952, 2383925, 2274286, 2234363, 2298409, 1722082, 2027459, 1959843, 2472748, 2229331, 2307136, 2040152, 2295857, 1943462, 1728326, 2050494, 2356113, 2064147, 1959875, 2110269, 1911993, 2481765, 1821588, 1905964, 2003304, 2100006, 2141153, 2028223, 2148632, 2237409, 2005323, 2021081, 2493279, 1992594, 1992514, 2308342, 2146914, 1776717, 2000576, 1839751, 2381455, 2198437, 2456263, 1867878, 2290853, 2385516]

selectiondec = np.where((starcat['_DEJ2000'] >= -90.0) & (starcat['_DEJ2000'] <= -30.0))
gaiaparallax = np.where(gaiacat['Parallax'][selectiondec] > 0)[0] #burp.gaiacat?
#absmag = burp.starcat['Vmag'][gaiaparallax]-5*np.log10(1/(gaiacat['parallax'][gaiaparallax]/1000)) + 5
absmag = starcat['Vmag'][gaiaparallax]-5*np.log10(1/(gaiacat['parallax'][gaiaparallax]/1000)) + 5
gaiamag = (starcat['Bmag']-starcat['Vmag'])[gaiaparallax]
#gaiamag = (burp.starcat['Bmag']-burp.starcat['Vmag'])[gaiaparallax]
print(gaiamag.size)
lim = np.where((gaiamag < 3) & (absmag < 10) & (absmag > -8))[0]
print(gaiamag[lim].size)


fig, ax = plt.subplots(1,1,figsize=(6,6))

h = ax.hist2d(gaiamag[lim], absmag[lim], bins=300, cmin=10, norm=colors.PowerNorm(0.5), zorder=0.5)
# fill the rest with scatter (set rasterized=True if saving as vector graphics)
ax.scatter(gaiamag[lim], absmag[lim], alpha=0.05, s=1, color='k', zorder=0, rasterized=True)
#h = ax.hist2d(gaiacat.bprp, absmag[lim], bins=300, cmin=10, norm=colors.PowerNorm(0.5), zorder=0.5)
# fill the rest with scatter (set rasterized=True if saving as vector graphics)
ax.scatter(gaiamag[lim], absmag[lim], alpha=0.05, s=1, color='k', zorder=0, rasterized=True)
ax.invert_yaxis()


#ax.plot(gaiamag[lim], absmag[lim],'.',color="#999999")
ax.plot(0.63, 4.83, 'o', color = 'gold')
ax.plot(0, 0.582, 'o', color = 'cyan')

ax.set_title("HR Diagram")
ax.set_xlabel("B-V Colour")
ax.set_ylabel("Absolute Magnitude")
fig.legend(['all bRing stars', 'Sun', 'Vega', 'mag < 4'], loc = "lower right")


plt.draw()
plt.savefig(paths.figures / 'hrd.pdf')
#plt.show()
print(starcat['ASCC'])
for i in stars:
    print(i)
    print(np.array(starcat['ASCC']) == i)
    quit()
    starcat2loc = np.where(starcat['ASCC'] == str(i))
#    starcat2loc = np.where(starcat['ASCC'] == i)[0][0]
    starinfo2 = starcat[starcat2loc]
#    starcat2loc = np.where(starcat2['ASCC'] == i)[0][0]
#    starinfo2 = starcat2[starcat2loc]
    par = starinfo2['Plx']
    vmag = starinfo2['Vmag']
    bmag = starinfo2['Bmag']
    absmag = vmag-5*np.log10(1/(par/1000)) + 5
    plt.plot(bmag - vmag, absmag, '.', color = 'r')

plt.gca().invert_yaxis()
plt.show()
