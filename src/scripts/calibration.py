import numpy as np
import sys
import astropy.stats as aps
import astropy.io.fits as fits
import astropy
import fromdaily
import fromquarterly
import background
import h5py
#from mpcubed.calibration import grids
import paths

#sys.path.append('/disks/web1/bring/products/')
#sys.path.append('/data2/stuik/code')

def openfits(filename):
  image_data = fits.getdata(filename)
  return image_data

starcat    = fits.getdata(paths.data / 'bringcat20191231.fits.gz')

class SysFile(object):
    # Read data from systematics files.
    #Attributes: sysfile (str): The full path to the file.
    
    def __init__(self, sysfile):
        # Initialize a reader of systematics files.
        #Args: sysfile (str): The full path to the file.
        
        self.sysfile = sysfile
        return

    def read_header(self):
        # Read all header attributes.
        #Returns: data: A list of attribute (key, value) pairs.
        
        with h5py.File(self.sysfile, 'r') as f:
            data = f['header'].attrs.items()
        return data

    def read_pointing(self):
        # Read the pointing associated with the systematics.
        #Returns:
        #    alt0 (float): The pointing altitude in degrees.
        #    az0 (float): The pointing azimuth in degrees.
        #    th0 (float): The pointing orientation in degrees.
        #    x0 (float): The pointing center in pixels.
        #    y0 (float): The pointing center in pixels.
        
        with h5py.File(self.sysfile, 'r') as f:
            alt0 = f['header'].attrs['alt0']
            az0 = f['header'].attrs['az0']
            th0 = f['header'].attrs['th0']
            x0 = f['header'].attrs['x0']
            y0 = f['header'].attrs['y0']
        return alt0, az0, th0, x0, y0

    def read_spatial(self):
        # Read the spatial quality statistics.
        #Returns: spatial (dict): A dictionary containing the quality statistics.
        
        spatial = dict()
        with h5py.File(self.sysfile, 'r') as f:
            grp = f['header/spatial']
            spatial['niter'] = grp['niter'][()]
            spatial['chisq'] = grp['chisq'][()]
            spatial['npoints'] = grp['npoints'][()]
            spatial['npars'] = grp['npars'][()]
        return spatial

    def read_temporal(self):
        # Read the temporal quality statistics.
        #Returns: temporal (dict): A dictionary containing the quality statistics.
        
        temporal = dict()
        with h5py.File(self.sysfile, 'r') as f:
            grp = f['header/temporal']
            temporal['niter'] = grp['niter'][()]
            temporal['chisq'] = grp['chisq'][()]
            temporal['npoints'] = grp['npoints'][()]
            temporal['npars'] = grp['npars'][()]
        return temporal

    def read_magnitudes(self):
        # Read the best-fit magnitudes.
       #Returns: magnitudes (dict): Containing: ascc, vmag, nobs, mag and sigma.
        
        magnitudes = dict()
        with h5py.File(self.sysfile, 'r') as f:
            grp = f['data/magnitudes']
            magnitudes['ascc'] = grp['ascc'][()]
            magnitudes['vmag'] = grp['vmag'][()]
            magnitudes['nobs'] = grp['nobs'][()]
            magnitudes['mag'] = grp['mag'][()]
            magnitudes['sigma'] = grp['sigma'][()]
        return magnitudes

    def read_trans(self):
        # Read the fitted transmission.
        #Returns:
        #    grid: A polar grid instance corresponding to the transmission map.
        #    trans (dict): Containing: grid, num_k, num_n, nobs, trans.
        
        trans = dict()
        with h5py.File(self.sysfile, 'r') as f:
            grp = f['data/trans']
            trans['gridtype'] = grp.attrs['grid']
            trans['num_k'] = grp.attrs['nx']
            trans['num_n'] = grp.attrs['ny']
            grid = grids.PolarGrid(trans['num_k'], trans['num_n'])
            idx1 = grp['idx1'][()]
            idx2 = grp['idx2'][()]
            trans['nobs'] = grid.values2grid(idx1, idx2, grp['nobs'][()])
            trans['trans'] = grid.values2grid(idx1, idx2, grp['trans'][()])
        return grid, trans

    def read_intrapix(self):
        # Read the fitted intrapixel variations.
        #Returns:
        #    grid: A polar grid instance corresponding to the intrapixel
        #        variations.
        #    intrapix (dict): Containing grid, num_l, num_n, nobs, amplitudes/
        
        intrapix = dict()
        with h5py.File(self.sysfile, 'r') as f:
            grp = f['data/intrapix']
            intrapix['gridtype'] = grp.attrs['grid']
            intrapix['num_l'] = grp.attrs['nx']
            intrapix['num_n'] = grp.attrs['ny']
            grid = grids.PolarGrid(intrapix['num_l'], intrapix['num_n'])
            idx1 = grp['idx1'][()]
            idx2 = grp['idx2'][()]
            sinx = grid.values2grid(idx1, idx2, grp['sinx'][()])
            cosx = grid.values2grid(idx1, idx2, grp['cosx'][()])
            siny = grid.values2grid(idx1, idx2, grp['siny'][()])
            cosy = grid.values2grid(idx1, idx2, grp['cosy'][()])
            intrapix['nobs'] = grid.values2grid(idx1, idx2, grp['nobs'][()])
            intrapix['amplitudes'] = np.stack([sinx, cosx, siny, cosy], axis=-1)
        return grid, intrapix

    def read_clouds(self):
        # Read the fitted clouds.
        #Returns:
        #    grid: A healpix grid instance corresponding to the clouds.
        #    clouds (dict): Containing: gridtype, num_q, lstmin, lstmax, lstlen,
        #        nobs, clouds, sigma.
        
        clouds = dict()
        with h5py.File(self.sysfile, 'r') as f:
            grp = f['data/clouds']
            clouds['gridtype'] = grp.attrs['grid']
            clouds['num_q'] = grp.attrs['nx']
            clouds['lstmin'] = grp.attrs['lstmin']
            clouds['lstmax'] = grp.attrs['lstmax']
            clouds['lstlen'] = grp.attrs['lstlen']
            grid = grids.HealpixGrid(clouds['num_q'])
            idx = grp['idx'][()]
            lstseq = grp['lstseq'][()] - clouds['lstmin']
            nobs_ = grp['nobs'][()]
            clouds_ = grp['clouds'][()]
            sigma_ = grp['sigma'][()]
        clouds['nobs'] = np.full((grid.npix, clouds['lstlen']), fill_value=np.nan)
        clouds['nobs'][idx, lstseq] = nobs_
        clouds['clouds'] = np.full((grid.npix, clouds['lstlen']), fill_value=np.nan)
        clouds['clouds'][idx, lstseq] = clouds_
        clouds['sigma'] = np.full((grid.npix, clouds['lstlen']), fill_value=np.nan)
        clouds['sigma'][idx, lstseq] = sigma_
        return grid, clouds
        
def get_calibration(ascc, dataD, dataQ, camera):
    basepath = '/data2/bring/reduced/' #was originally /data3/bring/reduced but cannot be found
    biweekly = ['A', 'B']
    indexfiles = np.uint32(np.unique(dataQ['findex']))
    datafiles = np.uint32(indexfiles/10)
    quarternames = np.array([str(np.uint32(ddd/100))+'Q'+str(np.uint32((ddd%100-1)/3+1)) for ddd in datafiles])
    filenames = np.array([str(np.uint32(iii/10))+biweekly[np.uint32((iii % 10)/5)] for iii in indexfiles])
    sysfiles = [basepath+qqq+'/'+camera+'/sys0_vmag_'+fff+camera+'.hdf5' for qqq,fff in zip(quarternames, filenames)]

    asccsel = np.where(starcat['ascc'] == str(ascc))[0][0]
    ra = starcat['_RAJ2000'][asccsel]
    dec = starcat['_DEJ2000'][asccsel]
    vmag = starcat['Vmag'][asccsel]
    incl = starcat['Inclusion'][asccsel]
    #ha = np.mod(ra - dataD['lst']*15., 360.)  # -- which one is right????
    ha = np.mod(dataD['lst']*15 - ra, 360.)
    decs = np.ones(len(ha))*dec
    ras = np.ones(len(ha))*ra
    allsys = np.zeros(len(ha), dtype='float32')
    alltrans = np.zeros(len(ha), dtype='float32')+np.nan
    nobstrans = np.zeros(len(ha), dtype='float32')
    allipx = np.zeros(len(ha), dtype='float32')+np.nan
    nobsipx = np.zeros(len(ha), dtype='float32')
    allclouds = np.zeros(len(ha), dtype='float32')+np.nan
    sigmaclouds = np.zeros(len(ha), dtype='float32')+np.nan
    nobsclouds = np.zeros(len(ha), dtype='float32')
    allsystematics = np.zeros(len(ha), dtype='float32')+np.nan
    cflag = np.zeros(len(ha), dtype='uint8')+32

    flux, eflux = dataD['flux0']/dataD['exptime'], dataD['eflux0']/dataD['exptime']
    sky, esky = dataD['sky']/dataD['exptime'], dataD['esky']/dataD['exptime']

    for sss, iii in zip(sysfiles, indexfiles):
        print("Running "+str(iii))
        sysfile = SysFile(sss)
        Qsel = np.where(dataQ['findex'] == iii)[0]
        Dsel = np.where((dataD['lstseq'] >= np.min(dataQ['lstseq'][Qsel]))  & (dataD['lstseq'] <= np.max(dataQ['lstseq'][Qsel])))[0]
        cflag[Dsel] = 0

        syscamgrid, systrans = sysfile.read_trans()
        sysipxgrid, sysintrapix = sysfile.read_intrapix()
        sysskygrid, sysclouds = sysfile.read_clouds()

        # Transmission
        k, n = syscamgrid.radec2idx(ha[Dsel], decs[Dsel])
        alltrans[Dsel] = systrans['trans'][k, n]
        nobstrans[Dsel] = systrans['nobs'][k, n]
        cflag[Dsel] = np.where(nobstrans[Dsel] < 25, cflag[Dsel]+2, cflag[Dsel])

        # Interpixel
        sinx, cosx = np.sin(2*np.pi*dataD['x'][Dsel]), np.cos(2*np.pi*dataD['x'][Dsel])
        siny, cosy = np.sin(2*np.pi*dataD['y'][Dsel]), np.cos(2*np.pi*dataD['y'][Dsel])
        b_mat = np.column_stack([sinx, cosx, siny, cosy])
        l, n = sysipxgrid.radec2idx(ha[Dsel], decs[Dsel])
        allipx[Dsel] = np.sum(b_mat*sysintrapix['amplitudes'][l,n], axis=1)
        nobsipx[Dsel] = sysintrapix['nobs'][l,n]
        cflag[Dsel] = np.where(nobsipx[Dsel] < 25, cflag[Dsel]+4, cflag[Dsel])

        # Clouds
        q = sysskygrid.radec2idx(ra, dec)
        t = dataD['lstseq'][Dsel] - sysclouds['lstmin']
        allclouds[Dsel] = sysclouds['clouds'][q,t]
        sigmaclouds[Dsel] = sysclouds['sigma'][q,t]
        nobsclouds[Dsel] = sysclouds['nobs'][q,t]
        cflag[Dsel] = np.where(nobsclouds[Dsel] < 25, cflag[Dsel]+8, cflag[Dsel])
        cflag[Dsel] = np.where(sigmaclouds[Dsel] > .05, cflag[Dsel]+16, cflag[Dsel])

        # Summing it all up
        allsystematics[Dsel] = alltrans[Dsel] + allipx[Dsel] + allclouds[Dsel]
        cflag[Dsel] = np.where(np.isnan(allsystematics[Dsel]), cflag[Dsel] + 1, cflag[Dsel])

    mag, emag = 25 - 2.5*np.log10(flux), 2.5/np.log(10.)*np.abs(eflux/flux)
    netmag = mag - alltrans - allipx - allclouds

    if len(np.where(np.array(dataD.dtype.names) == 'nmag0')[0]) == 0:
        dataD = np.lib.recfunctions.append_fields(dataD, 'nmag0', netmag, usemask = False)
    else:
        dataD['nmag0'] = netmag
    if len(np.where(np.array(dataD.dtype.names) == 'trans0')[0]) == 0:
        dataD = np.lib.recfunctions.append_fields(dataD, 'trans0', alltrans, usemask = False)
    else:
        dataD['trans0'] = alltrans
    if len(np.where(np.array(dataD.dtype.names) == 'ipx0')[0]) == 0:
        dataD = np.lib.recfunctions.append_fields(dataD, 'ipx0', allipx, usemask = False)
    else:
        dataD['ipx0'] = allipx
    if len(np.where(np.array(dataD.dtype.names) == 'clouds0')[0]) == 0:
        dataD = np.lib.recfunctions.append_fields(dataD, 'clouds0', allclouds, usemask = False)
    else:
        dataD['clouds0'] = allclouds
    if len(np.where(np.array(dataD.dtype.names) == 'eclouds0')[0]) == 0:
        dataD = np.lib.recfunctions.append_fields(dataD, 'eclouds0', sigmaclouds, usemask = False)
    else:
        dataD['eclouds0'] = sigmaclouds
    if len(np.where(np.array(dataD.dtype.names) == 'ccflag')[0]) == 0:
        dataD = np.lib.recfunctions.append_fields(dataD, 'ccflag', cflag, usemask = False)
    else:
        dataD['ccflag'] = cflag
    if len(np.where(np.array(dataD.dtype.names) == 'ntrans0')[0]) == 0:
        dataD = np.lib.recfunctions.append_fields(dataD, 'ntrans0', nobstrans, usemask = False)
    else:
        dataD['ntrans0'] = nobstrans
    if len(np.where(np.array(dataD.dtype.names) == 'nclouds0')[0]) == 0:
        dataD = np.lib.recfunctions.append_fields(dataD, 'nclouds0', nobsclouds, usemask = False)
    else:
        dataD['nclouds0'] = nobsclouds
    if len(np.where(np.array(dataD.dtype.names) == 'nipx0')[0]) == 0:
        dataD = np.lib.recfunctions.append_fields(dataD, 'nipx0', nobsipx, usemask = False)
    else:
        dataD['nipx0'] = nobsipx
    if len(np.where(np.array(dataD.dtype.names) == 'tmag0')[0]) == 0:
        dataD = np.lib.recfunctions.append_fields(dataD, 'tmag0', dataD['nmag0']+dataD['ipx0'])
    else:
        dataD['tmag0'] = dataD['nmag0']+dataD['ipx0']
    return dataD   

## Calibrating data of the star using all definitions and sunalt < -18
## Output: good daily data of star
def calibrate(star, camera = 'AUW'): #without 5 sigma mag outliers removed
    starcatloc = np.where(starcat['ascc'] == str(star))[0][0]
    starinfo = starcat[starcatloc]
    HD = starinfo['HD']  
    print(HD)
    caldaily = openfits(paths.data / f'{HD}_{camera}.fits')
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
    return datafinal
