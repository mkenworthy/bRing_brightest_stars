#!/usr/bin/python
#tr -d '\r' < getdataALL.py > tmp.py ; mv -f tmp.py getdataALL.py ; chmod 755 getdataALL.py

import numpy as np
import warnings
import glob
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py


def fromquarterly(inascc, camera):
    """
    fromquarterly() collects the data from star {inascc} of camera {camera}.
    This data is stored on /net/mascara0/data3/mascara/reduced/. The data consists
    of the lstseq (??), jd (julian date), nobs, findex, lst, exptime (exposure time),
    x (x position in pixels on the image), y (y position in pixels on the image),
    mag0 (V magnitude of the star), emag0 (error on the V magnitude of the star),
    sky, esky, trans0, etrans0, clouds0, eclouds0 and the camera.
    INPUT:
        inascc: string with the ascc name of the star
        camera: string with the abreviation of the camera
    OUTPUT:
        curves: recarray containing all elements mentioned above
    """
    # setting up empty lists to be filled later
    camstr = ['AUE', 'AUW', 'LSC', 'LSE', 'LSN', 'LSS', 'LSW', 'SAE', 'SAW', 'LPC', 'LPE', 'LPN', 'LPS', 'LPW']
    basepath = '/net/mascara0/data3/mascara/reduced/'
    lstseq = []
    jd = []
    nobs = []
    findex = []
    lst = []
    exptime = []
    x = []
    y = []
    mag0 = []  # Vmag
    emag0 = []
    sky = []
    esky = []
    trans0 = []
    etrans0 = []
    clouds0 = []
    eclouds0 = []
    cameras = []

    # finding the file and extracting the different properties wanted
    files = np.sort(glob.glob(basepath + '20*/' + camera + '/' + 'red0_vmag_??????????.hdf5'))
    for fff in files:
        curcam = np.sum([iii * (fff.find(camstr[iii]) != -1) for iii in range(len(camstr))])
        try:

            f = h5py.File(fff, 'r')
            index = fff[-15:-8]

            if index[-1] == 'A':
                fileindex = np.uint(index[:-1]) * np.uint(10)
            else:
                fileindex = np.uint(index[:-1]) * np.uint(10) + np.uint(5)
            try:
                if curcam <= 8:
                    stars = np.array(f['stars/ascc'])
                    if len(np.where(stars == bytes(inascc, 'utf-8'))[0] > 0):  # changed in py3
                        lc = f['lightcurves/' + str(inascc)][()].copy()
                        lstseq = np.append(lstseq, lc['lstseq'] * 50)
                        jd = np.append(jd, lc['jd'])
                        lst = np.append(lst, lc['lst'])
                        findex = np.append(findex, np.zeros(len(lc['jd']), dtype=np.uint) + fileindex)
                        exptime = np.append(exptime, lc['exptime'])
                        nobs = np.append(nobs, lc['nobs'])
                        x = np.append(x, lc['x'])
                        y = np.append(y, lc['y'])
                        mag0 = np.append(mag0, lc['mag0'])
                        emag0 = np.append(emag0, lc['emag0'])
                        sky = np.append(sky, lc['sky'])
                        esky = np.append(esky, lc['esky'])
                        trans0 = np.append(trans0, lc['trans0'])
                        etrans0 = np.append(etrans0, lc['etrans0'])
                        clouds0 = np.append(clouds0, lc['clouds0'])
                        eclouds0 = np.append(eclouds0, lc['eclouds0'])
                        cameras = np.append(cameras, 0.0 * lc['eclouds0'] + curcam)
                else:
                    stars = np.array(f['header/ascc'])
                    if len(np.where(stars == bytes(inascc, 'utf-8'))[0] > 0):  # changed in py3
                        lc = f['data/' + str(inascc)][()].copy()
                        lstseq = np.append(lstseq, lc['lstseq'] * 50)
                        jd = np.append(jd, lc['jdmid'])
                        lst = np.append(lst, lc['lst'])
                        findex = np.append(findex, np.zeros(len(lc['jdmid']), dtype=np.uint) + fileindex)
                        exptime = np.append(exptime, np.zeros(len(lc['jdmid']), dtype=np.float) + 6.4)
                        nobs = np.append(nobs, lc['nobs'])
                        x = np.append(x, lc['x'])
                        y = np.append(y, lc['y'])
                        mag0 = np.append(mag0, lc['mag0'])
                        emag0 = np.append(emag0, lc['emag0'])
                        sky = np.append(sky, lc['sky'])
                        esky = np.append(esky, lc['esky'])
                        trans0 = np.append(trans0, lc['trans0'])
                        etrans0 = np.append(etrans0, lc['etrans0'])
                        clouds0 = np.append(clouds0, lc['clouds0'])
                        eclouds0 = np.append(eclouds0, lc['eclouds0'])
                        cameras = np.append(cameras, 0.0 * lc['eclouds0'] + curcam)

                f.close()
            except KeyError:
                f.close()
        except IOError:
            {}
    # construct a recarray with all the properties
    names = ['lstseq', 'jd', 'lst', 'exptime', 'nobs', 'x', 'y', 'mag0', 'emag0', 'sky', 'esky', 'trans0', 'etrans0',
             'clouds0', 'eclouds0', 'findex', 'camera']
    formats = ['uint32', 'float64', 'float64', 'float32', 'uint64', 'float32', 'float32', 'float32', 'float32',
               'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'uint32']
    ndata = len(lstseq)
    curves = np.recarray(ndata, names=names, formats=formats)
    curves['lstseq'] = np.array(lstseq)
    curves['jd'] = np.array(jd)
    curves['lst'] = np.array(lst)
    curves['exptime'] = np.array(exptime)
    curves['nobs'] = np.array(nobs)
    curves['x'] = np.array(x)
    curves['y'] = np.array(y)
    curves['mag0'] = np.array(mag0)
    curves['emag0'] = np.array(emag0)
    curves['sky'] = np.array(sky)
    curves['esky'] = np.array(esky)
    curves['trans0'] = np.array(trans0)
    curves['etrans0'] = np.array(etrans0)
    curves['clouds0'] = np.array(clouds0)
    curves['eclouds0'] = np.array(eclouds0)
    curves['findex'] = np.array(findex)
    curves['camera'] = np.array(cameras)
    return (curves)


