#!/usr/bin/python
#tr -d '\r' < getdata.py > tmp.py ; mv -f tmp.py getdata.py ; chmod 755 getdata.py

import os,sys
import numpy as np
import h5py
import glob
import astropy.io.fits as pf

def fromquarterly(inascc, camera):
    camstr = ['AUE', 'AUW', 'LSC', 'LSE', 'LSN', 'LSS', 'LSW', 'SAE', 'SAW', 'LPC', 'LPE', 'LPN', 'LPS', 'LPW']
    basepath = '/net/belterwijde/data2/bring/reduced/'
    lstseq = []
    jd = []
    nobs = []
    findex = []
    lst = []
    exptime =[]
    x = []
    y = []
    mag0 = []
    emag0 = []
    sky = []
    esky = []
    trans0 = []
    etrans0 = []
    clouds0 = []
    eclouds0 = []
    starfrac = []
    dr = []
    moonalt = []
    moonph = []
    sunalt = []
    cameras = []

    files = np.sort(glob.glob(basepath+'*/'+camera+'/'+'red0*.hdf5'))
    for fff in files:
        curcam = np.sum([iii*(fff.find(camstr[iii]) != -1) for iii in range(len(camstr))])
        try:
            astroexists = False
            f = h5py.File(fff,'r')
            g = h5py.File(fff[:-21]+'fast_'+fff[-16:], 'r')
            index = fff[-15:-8]
            if index[-1] == 'A':
                fileindex = np.uint(index[:-1])*np.uint(10)
            else:
                fileindex = np.uint(index[:-1])*np.uint(10)+np.uint(5)
            try:
                key = g['astrometry']
                astroexists = True
            except KeyError:
                {}
            try:
                stars = np.array(f['stars/ascc'])
                if len(np.where(stars == bytes(inascc, 'utf-8'))[0] > 0):  # changed in py3
                    lc      = f['lightcurves/'+str(inascc)][()].copy()
                    lstseq  = np.append(lstseq,lc['lstseq']*50)
                    curlstseq = lc['lstseq']*50
                    jd      = np.append(jd,lc['jd'])
                    lst     = np.append(lst,lc['lst'])
                    findex  = np.append(findex,np.zeros(len(lc['jd']),dtype=np.uint)+fileindex)
                    exptime = np.append(exptime,lc['exptime'])
                    nobs    = np.append(nobs,lc['nobs'])
                    x       = np.append(x,lc['x'])
                    y       = np.append(y, lc['y'])
                    mag0    = np.append(mag0, lc['mag0'])
                    emag0   = np.append(emag0, lc['emag0'])
                    sky     = np.append(sky, lc['sky'])
                    esky    = np.append(esky, lc['esky'])
                    trans0   = np.append(trans0, lc['trans0'])
                    etrans0  = np.append(etrans0, lc['etrans0'])
                    clouds0  = np.append(clouds0, lc['clouds0'])
                    eclouds0 = np.append(eclouds0, lc['eclouds0'])
                    if astroexists:
                        index = np.searchsorted(g['astrometry/lstseq'][()], curlstseq)
                        dr       = np.append(dr, g['astrometry/dr'][()][index].copy())
                        starfrac = np.append(starfrac, g['astrometry/ustars'][()][index]/np.float32(g['astrometry/fstars'][()][index]))
                    else:
                        dr       = np.append(dr,0.0*lc['eclouds0'])
                        starfrac = np.append(starfrac,0.0*lc['eclouds0']+1.0)
                    index = np.searchsorted(g['station/lstseq'][()], curlstseq)
                    moonalt       = np.append(moonalt, g['station/moonalt'][()][index].copy())
                    moonph        = np.append(moonph, g['station/moonph'][()][index].copy())
                    sunalt        = np.append(sunalt, g['station/sunalt'][()][index].copy())
                    cameras = np.append(cameras, 0.0*lc['eclouds0']+curcam)
                f.close()
                g.close()
            except KeyError:
                f.close()
                g.close()
        except IOError:
            {}

    names   = ['lstseq',      'jd',     'lst', 'exptime',   'nobs',       'x',       'y',    'mag0',   'emag0',     'sky',    'esky',  'trans0', 'etrans0', 'clouds0', 'eclouds0',  'findex',      'dr',  'sunalt', 'moonalt', 'moonph', 'starfrac', 'camera']
    formats = ['uint32', 'float64', 'float64', 'float32', 'uint64', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32',  'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'uint32']
    ndata = len(lstseq)
    curves = np.recarray(ndata, names=names, formats=formats)
    curves['lstseq']   = np.array(lstseq)
    curves['jd']       = np.array(jd)
    curves['lst']      = np.array(lst)
    curves['exptime']  = np.array(exptime)
    curves['nobs']     = np.array(nobs)
    curves['x']        = np.array(x)
    curves['y']        = np.array(y)
    curves['mag0']     = np.array(mag0)
    curves['emag0']    = np.array(emag0)
    curves['sky']      = np.array(sky)
    curves['esky']     = np.array(esky)
    curves['trans0']   = np.array(trans0)
    curves['etrans0']  = np.array(etrans0)
    curves['clouds0']  = np.array(clouds0)
    curves['eclouds0'] = np.array(eclouds0)
    curves['findex']   = np.array(findex)
    curves['dr']       = np.array(dr)
    curves['moonalt']  = np.array(moonalt)
    curves['moonph']   = np.array(moonph)
    curves['sunalt']   = np.array(sunalt)
    curves['starfrac'] = np.array(starfrac)
    curves['camera']   = np.array(cameras)
    return curves

def fromdaily(inascc, camera, debug=False, getclouds=True):
    if debug:
        print('Looking for data on star: '+str(inascc) + ' in camera '+str(camera))
    camstr = ['AUE', 'AUW', 'LSC', 'LSE', 'LSN', 'LSS', 'LSW', 'SAE', 'SAW', 'LPC', 'LPE', 'LPN', 'LPS', 'LPW']
    # Get main data from light curves, combine with JD from station
    basepath = '/net/belterwijde/data2/bring/products/'
    if getclouds:
        data = pf.getdata('/net/belterwijde/data2/bring/Catalogues/bringcat20191231.fits')
        ascc = data['ASCC']
        skyidx = data['SKYIDX']
        sort = np.argsort(ascc)
        ascc = ascc[sort]
        skyidx = skyidx[sort]
        args = np.searchsorted(ascc, inascc)
        curskyidx = skyidx[args]
        if debug:
            print('Equivalent cloud IDX: '+str(curskyidx))

    lstseq = []
    jd = []
    lst = []
    exptime =[]
    x = []
    y = []
    flux0 = []
    flux1 = []
    eflux0 = []
    eflux1 = []
    mag0 = []
    emag0 = []
    cmag0 = []
    newsky = []
    peak = []
    sky = []
    esky = []
    aflag = []
    pflag = []
    cflag = []
    cameras = []
    clouds0 = []
    eclouds0 = []

    files = np.sort(glob.glob(basepath+'*'+camera+'/lightcurves/fast*'))
    if debug:
        print('Found '+str(len(files))+' datafiles meeting criteria.')
    for fff in files:
        curcam = np.sum([iii*(fff.find(camstr[iii]) != -1) for iii in range(len(camstr))])
        try:
            f = h5py.File(fff, 'r')
            try:
                stars = np.array(f['stars/ascc'])
                if len(np.where(stars == bytes(inascc, 'utf-8'))[0] > 0):  # changed in py3
                    if debug:
                        print('Found data in ' + fff)
                    station = f['station']
                    if len(np.where([iii == 'LSTSEQ' for iii in station.keys()])[0] > 0):  # changed in py3
                        inlstseq  = station['LSTSEQ'][()].copy()
                        injd      = station['JD'][()].copy()
                        inlst     = station['LST'][()].copy()
                        inexptime = station['EXPTIME'][()].copy()
                    else:
                        inlstseq  = station['lstseq'][()].copy()
                        injd      = station['jd'][()].copy()
                        inlst     = station['lst'][()].copy()
                        inexptime = station['exptime'][()].copy()
                    lc        = f['lightcurves/'+str(inascc)][()].copy()
                    curlst    = lc['lstseq']
                    idx       = np.searchsorted(inlstseq,curlst)

                    lstseq    = np.append(lstseq,curlst)
                    jd        = np.append(jd,injd[idx])
                    lst       = np.append(lst,inlst[idx])
                    exptime   = np.append(exptime,inexptime[idx])
                    x         = np.append(x,lc['x'])
                    y         = np.append(y,lc['y'])
                    flux0     = np.append(flux0,lc['flux0'])
                    flux1     = np.append(flux1,lc['flux1'])
                    eflux0    = np.append(eflux0,lc['eflux0'])
                    eflux1    = np.append(eflux1,lc['eflux1'])
                    peak      = np.append(peak,lc['peak'])
                    sky       = np.append(sky,lc['sky'])
                    esky      = np.append(esky,lc['esky'])
                    aflag     = np.append(aflag,lc['aflag'])
                    pflag     = np.append(pflag,lc['pflag'])
                    tempsky   = np.nan*np.zeros(len(lc['sky']))
                    tempcmag0 = np.zeros(len(lc['sky']))
                    cameras   = np.append(cameras, 0.0*lc['sky']+curcam)
                    if len(np.where([iii == 'tmpmag0' for iii in lc.dtype.names])[0] > 0):  # changed in py3
                        mag0    = np.append(mag0,lc['tmpmag0'])
                        emag0   = np.append(emag0,lc['tmpemag0'])
                        cflag   = np.append(cflag,lc['cflag'])
                    else:
                        mag0    = np.append(mag0,np.nan*np.zeros(len(lc['sky'])))
                        emag0   = np.append(emag0,np.nan*np.zeros(len(lc['sky'])))
                        cflag   = np.append(cflag,np.nan*np.zeros(len(lc['sky'])))
                    # Try to add the calibration magnitudes for each star
                    try:
                        split = fff.split('fast')
                        mcurlst = np.floor(curlst / 50)*50
                        corrfile = split[0]+'corr'+split[1]
                        g = h5py.File(corrfile, 'r')
                        stars = np.array(g['lightcurves'].keys())
                        if len(np.where(stars == inascc)[0] > 0):
                            lc      = g['lightcurves/'+str(inascc)][()].copy()
                            modlst = np.floor(lc['lstseq'] / 50)*50
                            for mmm, sss, ccc in zip(modlst, lc['sky'], lc['cmag0']):
                                setidx = np.where(mcurlst == mmm)[0]
                                if len(setidx) > 0:
                                    tempsky[setidx] = sss
                                    tempcmag0[setidx] = ccc
                        g.close()
                        if debug:
                            print('** Found corrected magnitude data')
                        newsky = np.append(newsky, tempsky)
                        cmag0 = np.append(cmag0, tempcmag0)
                    except IOError:
                        {}
                        newsky = np.append(newsky, np.nan*np.zeros(len(lc['sky'])))
                        cmag0 = np.append(cmag0, np.nan*np.zeros(len(lc['sky'])))
                        f.close()
                    if getclouds:
                        try:
                            split = fff.split('lightcurves')
                            sysfile = split[0]+'sys/sys0_vmag'+split[1][5:]
                            h = h5py.File(sysfile, 'r')
                            clouds = h['data/clouds']
                            select = np.where(np.array(clouds['idx']) == curskyidx)[0]
                            lstseqarr = np.array(clouds['lstseq'])
                            cloudsarr = np.array(clouds['clouds'])
                            ecloudsarr = np.array(clouds['sigma'])
                            tmpclouds = []
                            tmpeclouds = []
                            db = {}
                            for sss in select:
                                db[lstseqarr[sss]] = [cloudsarr[sss], ecloudsarr[sss]]
                            for lll in curlst:
                                if lll in db:
                                    tmpclouds.append(db[lll][0])
                                    tmpeclouds.append(db[lll][1])
                                else:
                                    tmpclouds.append(np.nan)
                                    tmpeclouds.append(np.nan)
                            h.close()
                            if debug:
                                print('** Found clouds data')
                            clouds0 = np.append(clouds0, tmpclouds)
                            eclouds0 = np.append(eclouds0, tmpeclouds)
                        except IOError:
                            {}
                            clouds0 = np.append(clouds0, np.nan*np.zeros(len(lc['sky'])))
                            eclouds0 = np.append(eclouds0, np.nan*np.zeros(len(lc['sky'])))
                            f.close()
            except KeyError:
                f.close()
        except IOError:
            {}
    if debug:
        print('Data mining done, combining into recarray.')

    if getclouds:
        names   = ['lstseq', 'jd', 'lst', 'exptime', 'x', 'y', 'flux0', 'flux1', 'eflux0', 'eflux1', 'mag0', 'emag0', 'peak', 'sky', 'esky', 'aflag', 'pflag', 'cflag', 'cmag0','newsky', 'camera', 'clouds0', 'eclouds0']
        formats = ['uint32', 'float64', 'float64', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'uint16', 'uint16', 'uint16', 'float32', 'float64', 'uint32', 'float64', 'float64']
    else:
        names   = ['lstseq', 'jd', 'lst', 'exptime', 'x', 'y', 'flux0', 'flux1', 'eflux0', 'eflux1', 'mag0', 'emag0', 'peak', 'sky', 'esky', 'aflag', 'pflag', 'cflag', 'cmag0','newsky', 'camera']
        formats = ['uint32', 'float64', 'float64', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'uint16', 'uint16', 'uint16', 'float32', 'float64', 'uint32']
    ndata = len(lstseq)
    curves = np.recarray(ndata, names=names, formats=formats)
    curves['lstseq']  = np.array(lstseq)
    curves['jd']      = np.array(jd)
    curves['lst']     = np.array(lst)
    curves['exptime'] = np.array(exptime)
    curves['x']       = np.array(x)
    curves['y']       = np.array(y)
    curves['flux0']   = np.array(flux0)
    curves['flux1']   = np.array(flux1)
    curves['eflux0']  = np.array(eflux0)
    curves['eflux1']  = np.array(eflux1)
    curves['mag0']    = np.array(mag0)
    curves['emag0']   = np.array(emag0)
    curves['cmag0']   = np.array(cmag0)
    curves['newsky']  = np.array(newsky)
    curves['peak']    = np.array(peak)
    curves['sky']     = np.array(sky)
    curves['esky']    = np.array(esky)
    curves['aflag']   = np.array(aflag)
    curves['pflag']   = np.array(pflag)
    curves['cflag']   = np.array(cflag)
    curves['camera']  = np.array(cameras)
    if getclouds:
        curves['clouds0']  = np.array(clouds0)
        curves['eclouds0']  = np.array(eclouds0)
    return curves
