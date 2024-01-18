import os
import numpy as np
#import pyfits
import astropy.io.fits as pyfits

def mrdfits(filename, extension=1, transpose=False, verbose=False):
    """
    Read a an extended fits file.
    Try to emulate mrdfits.pro

    Arguments
    ---------
    filename: string
    extension: int
    transpose: bool
    verbose: bool

    Returns
    -------
    A dictionary of arrays

    Update:
    ------
    fitsopen changed to open 20 january 2012 Clerc Gastaud

    Mer 11 jan 2012 data._names change to data.names


    """
    if verbose:
        print('my extension = ', extension)
    # load data
    fits = pyfits.open(filename)
    hdu = fits[extension]
    data = hdu.data
    if not (isinstance(hdu, pyfits.BinTableHDU)
            or isinstance(hdu, pyfits.TableHDU)):
        if verbose: print( 'this is an "image"')
        #print 'data', data.dtype, data.size
        if ((data.dtype == np.uint8) and (data.size > 1)):
            data = np.ndarray.tostring(data)
            print('conversion to string of data=', data)
        return data
    # else this is a table (ASCII or Binary)
    header = hdu.header
    #  create dictionnary
    mydict = dict.fromkeys(data.names)
    # beware order is not kept
    # fill the dictionnary
    for i, name in enumerate(data.names):
       mydim = "tdim{0}".format(i+1)
       #if header.has_key(mydim):
       if mydim in header:
           myshape = header[mydim]
           # n = [int(i) for i in header.strip('()').split(',')]
           exec('myshape = ' + header[mydim])
           # need to reverse shape due to FITS convention
           # (quick axis first in FITS and last in Python)
           myshape = myshape[::-1]
           mydict[name] = data.field(name).reshape(myshape)
           if transpose:
               mydict[name] = mydict[name].transpose()
       else:
           mydict[name] = data.field(name).ravel()
    #  15 May 2015
    #  strings are stored as byte arrays pb of pyfits
    for k, v in mydict.items():
        #print 'v=', v.dtype, v.size, v
        if ((v.dtype == np.uint8) and (v.size > 1)):
                mydict[k] = np.ndarray.tostring(v)
                print('conversion to string of v=', v)

    # check
    if verbose:
        for k, v in mydict.items():
            if  not isinstance(v, (str)):
                print(k, v[0:1], v.dtype)
            else:
                print(k, v)
    #
    return mydict

def mwrfits(mydict, filename, clobber=False, ascii=False, append=False, extname=None):
    """
    Write an dictionary as a binary table in a FITS file.

    Shape order is reverted to satisfy the FITS conventions
    (quick axis is first in FITS and second in python)

    Arguments
    ---------
    mydict: dict of ndarrays
    filename: string
      The name of the written fits file.
    clobber: bool
       Overwrites the output file if True.

    Returns
    --------
    Returns nothing.
    
    Exceptions
    ---------
    ValueError if mydict is not a dict of arrays
    """
    # check that dict contains only arrays
    for key, value in mydict.items():
        if not isinstance(value, np.ndarray):
            print('before key, type',key, type(value))
            if  not isinstance(value, (str)):
                value = np.asarray(value).reshape(1)
            else:
                value = np.fromstring(value, dtype=np.uint8)
            mydict[key] = value
            print('after key, type',key, type(value))
            print('')
            #raise ValueError("Expected a dict of arrays.")

    # convert dict of ndarray as an array of size 1
    mydtype = list()
    for k, v in mydict.items():
        mydtype.append((k, str(v.shape) + v.dtype.str))
    arr = np.zeros(1, dtype=mydtype)
    for k, v in mydict.items():
        arr[k] = v
    if ascii:
        raise NotImplemented() # hdu = pyfits.TableHDU(arr) # not working !
    else:
        hdu = pyfits.BinTableHDU(arr)

    if extname is not None: hdu.header['EXTNAME'] = extname
    
    # shape order is reverted to satisfy the FITS conventions
    # (quick axis is first in FITS and second in python)
    for i, v in enumerate(mydict.values()):
        hdu.header['TDIM' + str(i + 1)] = str(v.shape[::-1])
    if append:
        print(' append ')
        fits = pyfits.open(filename)
        fits.append(hdu)
        fits.writeto(filename, clobber=True)
    else:
        hdu.writeto(filename, clobber=clobber)
