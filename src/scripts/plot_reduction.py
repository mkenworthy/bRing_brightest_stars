#plot all steps of outlier rejection
#cloud error outliers were removed before we saved Remkos raw data
import astropy.io.fits as fits

import matplotlib.pyplot as plt
import numpy as np
import background

import paths

starcat    = fits.getdata(paths.data / 'bringcat20191231.fits.gz')

def openfits(filename):
  image_data = fits.getdata(filename)
  return image_data


#pick a star to show data reduction steps
star = 2424942

#calibration.calibrate
starcatloc = np.where(starcat['ascc'] == str(star))[0][0]
starinfo = starcat[starcatloc]
HD = starinfo['HD']  
camera = 'AUW' # RawData
#caldaily = openfits('RawData/'+str(HD)+'_'+camera+'.fits')
file = str(HD)+'_'+camera+'.fits'
print(file)
caldaily = openfits(paths.data / file)
if camera == 'AUE': #remove data with unknown magnitude outliers for camera AUE
    mask = np.where((caldaily["jd"] < 2458212) | (caldaily["jd"] >= 2458214))
    caldaily = caldaily[mask[0]]
elif camera == 'SAE': #remove data with unknown magnitude outliers for camera SAE
    mask = np.where((caldaily["jd"] < 2459103) | (caldaily["jd"] >= 2459104))
    caldaily = caldaily[mask[0]]
clouddaily = background.dailycloud(caldaily)
skydaily = background.sky(clouddaily)
magdaily = background.dailymag(skydaily)
locsun = (magdaily['sunalt'] <= -18) 
data = magdaily[locsun]
prepreselect = np.where( ((data['aflag'] <= 1) | (data['aflag'] == 0)) & ((data['ccflag'] == 0) | (data['ccflag'] == 4)) & (np.isfinite(data['tmag0'])) & (np.isfinite(data['emag0'])) & (data['emag0'] != 0) & (data['sunalt'] < -18) )[0]
datafinal = data[prepreselect]

fig, axs = plt.subplots(6, figsize=(12,12),sharex= True, sharey = False)

stu = dict(alpha=0.1, rasterized=True,marker='o',mec='none',markersize=1,color='brown',linestyle="")

mj=2400000.5
axs[0].plot(caldaily['jd']-mj, caldaily['mag0'], **stu)
axs[1].plot(clouddaily['jd']-mj, clouddaily['mag0'], **stu)
axs[2].plot(skydaily['jd']-mj, skydaily['mag0'], **stu)
axs[3].plot(magdaily['jd']-mj, magdaily['mag0'], **stu)
axs[4].plot(data['jd']-mj, data['mag0'], **stu)
axs[5].plot(datafinal['jd']-mj, datafinal['mag0'], **stu)
# axs[0].set(ylabel = "Magnitude")
# axs[1].set(ylabel = "Magnitude")
# axs[2].set(ylabel = "Magnitude")
# axs[3].set(ylabel = "Magnitude")
# axs[4].set(ylabel = "Magnitude")
fig.supylabel('Magnitude',fontsize=20)

axs[5].set_xlabel("Time (Modified Julian Date)",fontsize=20)
axs[0].set_title("Raw data, short exposures only, cloud error outliers rejected")
axs[1].set_title("Cloud outliers rejected")
axs[2].set_title("Sky background outliers rejected")
axs[3].set_title("Magnitude error outliers rejected")
axs[4].set_title("Sun altitude selection")
axs[5].set_title("Bad astrometry and photometry flags rejected, only finite magnitude values")
axs[0].invert_yaxis()
axs[1].invert_yaxis()
axs[2].invert_yaxis()
axs[3].invert_yaxis()
axs[4].invert_yaxis()
axs[5].invert_yaxis()
axs[5].set_ylim(3.8,3.3)
fig.suptitle('ASCC '+ str(star)+ ', HD ' + str(HD))
plt.draw()
plt.savefig(paths.figures / 'reduction.pdf')
#plt.show()


