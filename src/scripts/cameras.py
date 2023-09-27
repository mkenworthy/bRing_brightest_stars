import ephem
import numpy as np

jd1        = 2457906.0                      # June 1, 2017 -- date where sky was recalibrated
cams       = ['AUE', 'AUW', 'LSC', 'LSE', 'LSN', 'LSS', 'LSW', 'SAE', 'SAW', 'LPC', 'LPE', 'LPN', 'LPS', 'LPW', 'ESS', 'ESN']
obsSA = ephem.Observer()
obsSA.long = 20.8102*np.pi/180
obsSA.lat = -32.3812*np.pi/180
obsSA.elev = 1798.
obsSA.date = jd1-2415020
obsSA.epoch = obsSA.date

obsAU = ephem.Observer()
obsAU.long = 149.062150*np.pi/180
obsAU.lat = -31.272189*np.pi/180
obsAU.elev = 1165.
obsAU.date = jd1-2415020
obsAU.epoch = obsAU.date

obsLS = ephem.Observer()
obsLS.long = -70.727795*np.pi/180
obsLS.lat = -29.257124*np.pi/180
obsLS.elev = 2369.4
obsLS.date = jd1-2415020
obsLS.epoch = obsLS.date

obsLP = ephem.Observer()
obsLP.long = -17.8792*np.pi/180
obsLP.lat = 28.76025*np.pi/180
obsLP.elev = 2364.
obsLP.date = jd1-2415020
obsLP.epoch = obsLP.date

obsESS = ephem.Observer()
obsESS.long = -70.8059205*np.pi/180
obsESS.lat = -30.167622*np.pi/180
obsESS.elev = 2200.
obsESS.date = jd1-2415020
obsESS.epoch = obsLP.date

obsESN = ephem.Observer()
obsESN.long = -116.428161*np.pi/180
obsESN.lat = 32.840866*np.pi/180
obsESN.elev = 1859.
obsESN.date = jd1-2415020
obsESN.epoch = obsLP.date

locsPE = [obsAU, obsAU, obsLS, obsLS, obsLS, obsLS, obsLS, obsSA, obsSA, obsLP, obsLP, obsLP, obsLP, obsLP, obsESS, obsESN]





