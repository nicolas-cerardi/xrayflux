import os
#atomdb = os.environ['ATOMDB']
#work = os.environ['WORK']
import pyatomdb
import numpy as np
from tqdm import tqdm
from astropy.io import fits as F
from scipy.interpolate import interp1d, interp2d
from xrayflux.fitsutils import mwrfits
import matplotlib.pyplot as plt

def compute_flux(range_E=[0.5, 2.0], n_E=10, 
             T=1, 
             z=0., 
             abund=0.3,
             idx_start_abund=1,
             elements=None,
             rmf=None, 
             arf=None,
             eebrems=False,
             verbose=0,
             plots=False):
    
    ebins_in, ebins_out, ebins_in_m, ebins_out_m, rmfmatrix, arf, abres, elements = \
                                prepare_session_rmf_arf_elts(range_E, n_E, rmf, arf, elements, verbose=verbose)
    
    tmp_sess = pyatomdb.spectrum.CIESession(elements=elements)
    tmp_sess.set_eebrems(eebrems)
    if elements[0]!=1:
        # it's all metals anyway
        tmp_sess.set_abund(elements,abund)
        abres[:] = abund
    elif len(elements)>1:
        # whatever you ask
        tmp_sess.set_abund(elements[idx_start_abund:],abund)
        abres[idx_start_abund:] = abund
    # other case is just H alone, don't change its abundancy
    if plots:
        fig, axs = plt.subplots(1, 2, figsize=[15,10])
    tmp_ebins = ebins_in*(1+z) 
    tmp_sess.set_response(tmp_ebins, raw=True)

    # returns f/ne so in ph cm^3 s-1 bin-1
    tmp_spec = tmp_sess.return_spectrum(T)
    if plots: axs[0].loglog(mid_vals(ebins_in), tmp_spec, label=('%.1f keV'%T))
    # res in ph cm^3 s-1, or ph cm^5 s-1 if rmf and arf are used.
    if verbose>0 : print(tmp_spec.shape, rmfmatrix.shape)
    #apply rmf then arf
    if rmf is not None: tmp_spec = np.matmul(tmp_spec,rmfmatrix)
    tmp_spec *= arf #np.ones(ebins_out_m.shape) #
    if plots: axs[1].loglog(mid_vals(ebins_out), tmp_spec, label=('%.1f keV'%T))
    #integrate spectrum, apply l.o.s. dilatation and store it at the right (z,kT) index
    res = tmp_spec.sum()/(1+z) 
    if plots:
        axs[1].annotate('z=%.1f'%z, (0.2,0.2), xycoords='axes fraction', size=15, color='k')
        axs[1].legend()
        axs[0].set_xlim([range_E[0], range_E[1]])
        axs[1].set_xlabel('keV')
        axs[0].set_ylabel('ph cm^3 s^-1')
        axs[1].set_ylabel('ph cm^5 s^-1')
    return res
    
    
def prepare_session_rmf_arf_elts(range_E, n_E, rmf, arf, elements, verbose=0):
    # ebins in and out are because of the rmf. ebins in must be larger as some photons outside 
    # the ebins_out may be falsely detected inside ebins_out... we have to take account of them.
    # boundaries
    ebins_out = np.linspace(range_E[0],range_E[1],n_E)
    ebins_in = ebins_out.copy()
    # central values of bins
    ebins_out_m = mid_vals(ebins_out)
    ebins_in_m = mid_vals(ebins_in)
    if verbose>0:
        print(ebins_out)
        print(ebins_in)
    if rmf is not None:
        assert arf is not None
        ebins_in = np.linspace(range_E[0]*0.5,range_E[1]*2,n_E*4) #only useful in rmf is not None
        ebins_in_m = mid_vals(ebins_in)
        # load and prepare the aeff, arf, and their E range.
        rmfmatrix, arf = load_rmf_arf(rmf, arf, mid_vals(ebins_in), mid_vals(ebins_out))
    else : 
        # just a dummy aeff
        arf = np.ones(ebins_out_m.shape)
        rmfmatrix = None # this will make the code crash
    # abundancies and elements
    elements = elements if elements is not None else range(1, 3)
    abres = np.ones(np.array(elements).shape)
    return ebins_in, ebins_out, ebins_in_m, ebins_out_m, rmfmatrix, arf, abres, elements

def make_fluxtable(outfile=None,
             range_E=[0.5, 2.0], n_E=10, 
             range_T=[0.05, 50], n_T=10, 
             range_z=[0.0, 3.0], n_z=10, 
             abund=0.3,
             idx_start_abund=1,
             elements=None,
             rmf=None, 
             arf=None,
             eebrems=False,
             verbose=0,
             plots=False):
    """
    Generate tabulated xray flux data.

    Parameters:
    - outfile (str): Output file path.
    - range_E (list): Energy range.
    - n_E (int): Number of energy bins.
    - range_T (list): Temperature range.
    - n_T (int): Number of temperature bins.
    - range_z (list): Redshift range.
    - n_z (int): Number of redshift bins.
    - abund (float): Abundance value.
    - idx_start_abund (int): Index to start setting abundance.
    - elements (list): List of elements.
    - rmf (str): RMF file path.
    - arf (str): ARF file path.
    - eebrems (bool): Include electron-electron bremsstrahlung.
    - verbose (int): Verbosity level.
    - plots (bool): Generate plots.

    Returns:
    - tuple: Result, redshift bins, temperature bins, elements, abundance results.
    """
    # redshift and Temperature
    zbins = np.linspace(range_z[0],range_z[1],n_z)
    Tbins = np.logspace(np.log10(range_T[0]),np.log10(range_T[1]),n_T)
    
    if verbose>0:
        print(zbins)
        print(Tbins)
    res = np.zeros((n_z, n_T))
    
    ebins_in, ebins_out, ebins_in_m, ebins_out_m, rmfmatrix, arf, abres, elements = \
                                prepare_session_rmf_arf_elts(range_E, n_E, rmf, arf, elements, verbose=verbose)
    
    if plots:
        fig, axs = plt.subplots(2, n_z, figsize=[15,10])
        [[a.set_xscale('log') for a in ax] for ax in axs]
        [[a.set_yscale('log') for a in ax] for ax in axs]
        
    with tqdm(total=n_z * n_T) as pbar:
        tmp_sess = pyatomdb.spectrum.CIESession(elements=elements)
        tmp_sess.set_eebrems(eebrems)
        if elements[0]!=1:
            # it's all metals anyway
            tmp_sess.set_abund(elements,abund)
            abres[:] = abund
        elif len(elements)>1:
            # whatever you ask
            tmp_sess.set_abund(elements[idx_start_abund:],abund)
            abres[idx_start_abund:] = abund
        # other case is just H alone, don't change its abundancy
        
        for i, z in enumerate(zbins):
            # apply blueshift to ebins
            tmp_ebins = ebins_in*(1+z) 
            tmp_sess.set_response(tmp_ebins, raw=True)
            for j, kT in enumerate(Tbins):
                # returns f/ne so in ph cm^3 s-1 bin-1
                tmp_spec = tmp_sess.return_spectrum(kT)
                if plots: axs[0, i].plot(mid_vals(ebins_in), tmp_spec, label=('%.1f keV'%kT))
                # res in ph cm^3 s-1, or ph cm^5 s-1 if rmf and arf are used.
                if verbose>0 : print(tmp_spec.shape, rmfmatrix.shape)
                # apply rmf and arf
                if rmf is not None: tmp_spec = np.matmul(tmp_spec,rmfmatrix)
                tmp_spec *= arf
                if plots: axs[1, i].plot(mid_vals(ebins_out), tmp_spec/(1+z) , label=('%.1f keV'%kT))
                #integrate spectrum, apply l.o.s. dilatation and store it at the right (z,kT) index
                res[i,j] = tmp_spec.sum()/(1+z) 
                pbar.update(1)
            if plots:
                axs[0,i].annotate('z=%.1f'%z, (0.1,0.1), xycoords='axes fraction', size=15, color='k')
                axs[1,i].legend()
                axs[0,i].set_ylim([1e-20, None])
                axs[0,i].set_xlim(range_E)
                axs[1,i].set_xlim(range_E)
                axs[1,i].set_ylim([1e-20, None])
                axs[1,i].set_xlabel('keV')
        if plots:
            axs[0,0].set_ylabel(r'$ph\ cm^{3}\ s^{-1}\ bin^{-1}$')
            axs[1,0].set_ylabel(r'$ph\ cm^{5}\ s^{-1}\ bin^{-1}$')
    #print(elements, abres)
    if outfile is not None :
        outdict = {
            "TABLE":res,
            "RANGE_Z":zbins,
            "RANGE_T":Tbins,
            "ABUNDANCY":abres,
            "ELEMENTS":np.array(elements)
        }
        mwrfits(outdict,outfile,clobber=True)   
        
    return res, zbins, Tbins, elements, abres

def get_response_ebins(rmfdat):
    """
    mainly copied/adapted from https://github.com/AtomDB/pyatomdb/blob/c90494674dd01f3be204200238dac3fab1bd0af6/pyatomdb/pyatomdb/spectrum.py
    """
    matrixname='MATRIX'
    specbins_in = rmfdat[matrixname].data['ENERG_LO']
    specbins_in = np.append(specbins_in, rmfdat[matrixname].data['ENERG_HI'][-1])
    specbins_out = rmfdat['EBOUNDS'].data['E_MIN']
    specbins_out = np.append(specbins_out , rmfdat['EBOUNDS'].data['E_MAX'][-1])
    return specbins_in, specbins_out

def load_rmf_arf(rmf, arf,
                 specbins_in, specbins_out):
    """
    mainly copied/adapted from https://github.com/AtomDB/pyatomdb/blob/c90494674dd01f3be204200238dac3fab1bd0af6/pyatomdb/pyatomdb/spectrum.py
    
    returns an rmf matrix and a arf array. The loaded rmf is adjusted and interpolated from (specbins, ebins_out), 
    given by the rmf original file, to (specbins_in, specbins_out) asked by the code.
    Same for arf, interpolated into specbins_out.
    This clips energies wide outside the interesting range
    """
    rmf= F.open(rmf)
    arfdat = F.open(arf)
    arf = np.array(arfdat['SPECRESP'].data['SPECRESP'])
    ebins = rmf['EBOUNDS'].data['E_MIN']
    if ebins[-1] > ebins[0]:
        ebins = np.append(ebins, rmf['EBOUNDS'].data['E_MAX'][-1])
    else:
        ebins = np.append(rmf['EBOUNDS'].data['E_MAX'][0],ebins)
    matrixname = 'MATRIX'
    chanoffset = rmf['EBOUNDS'].data['CHANNEL'][0]

    rmfmatrix = np.zeros([len(rmf[matrixname].data),len(rmf['EBOUNDS'].data)])
    for ibin, i in enumerate(rmf[matrixname].data):
        lobound = 0
        fchan = i['F_CHAN']*1
        nchan = i['N_CHAN']*1
        if np.isscalar(fchan):
            fchan = np.array([fchan])
        fchan -= chanoffset
        if np.isscalar(nchan):
            nchan = np.array([nchan])
        for j in range(len(fchan)):
            ilo = fchan[j]
            if ilo < 0: continue
            ihi = fchan[j] + nchan[j]
            rmfmatrix[ibin,ilo:ihi]=i['MATRIX'][lobound:lobound+nchan[j]]
            lobound = lobound+nchan[j]
    specbins, ebins_out = get_response_ebins(rmf)

    if ebins_out[-1] < ebins_out[0]:
        # need to reverse things
        ebins_out=ebins_out[::-1]
        rmfmatrix = rmfmatrix[:,::-1]
    
    # find the boundaries of the interesting region of the rmf
    out_min, out_max = np.argwhere(ebins_out<specbins_out[0])[-1][0], \
                       np.argwhere(ebins_out>specbins_out[-1])[0][0]
    in_min, in_max = np.argwhere(specbins<specbins_in[0])[-1][0], \
                     np.argwhere(specbins>specbins_in[-1])[0][0]
    
    # Then perform the interpolation
    rmf_interpolated = interp2d(mid_vals(ebins_out[out_min:out_max+1]), mid_vals(specbins[in_min:in_max+1]), \
                               rmfmatrix[in_min:in_max,out_min:out_max])(specbins_out, specbins_in)
    
    rmf_interpolated = rmf_interpolated * mid_vals(ebins_out[out_min:out_max+1]).size/specbins_out.size
    
    arf_interpolated = interp1d(mid_vals(specbins), arf)(specbins_out)
    
    return rmf_interpolated, arf_interpolated

def mid_vals(arr):
    return (arr[1:]+arr[:-1])/2.0
