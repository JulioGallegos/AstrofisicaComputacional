def region_around_line(w, flux, cont):
    """cut out and normalize flux around a line
        
        Parameters
        ----------
        w : 1 dim np.ndarray
        array of wavelengths
        flux : np.ndarray of shape (N, len(w))
        array of flux values for different spectra in the series
        cont : list of lists
        wavelengths for continuum normalization [[low1,up1],[low2, up2]]
        that described two areas on both sides of the line
        """
    import numpy as np
    
    #index is true in the region where we fit the polynomial
    
    indcont = ((w > cont[0][0]) & (w < cont[0][1])) |((w > cont[1][0]) & (w < cont[1][1]))
    
    #index of the region we want to return
    
    indrange = (w > cont[0][0]) & (w < cont[1][1])
    
    # make a flux array of shape
    # (number of spectra, number of points in indrange)
    
    f = np.zeros((flux.shape[0], indrange.sum()))
    
    for i in range(flux.shape[0]):
        # fit polynomial of second order to the continuum region
        linecoeff = np.polyfit(w[indcont], flux[i, indcont],2)
        # divide the flux by the polynomial and put the result in our
        # new flux array
        f[i,:] = flux[i,indrange]/np.polyval(linecoeff, w[indrange])
    return w[indrange], f

def read_spec(filename):
    """ Read a UVES spectrum from the ESO pipeline
    
    Parameters
    ----------
    filename : string
    name of the fits file with the data
    
    Returns
    -------
    wavelength : np.ndarray
    wavelength (in Ang)
    flux : np.ndarray
    flux (in erg/s/cm**2)
    date_obs : string
    time of observation 
    """
    #Import modules
    
    from astropy.io import fits
    import numpy as np
    from astropy.wcs import WCS
    
    
    
    sp = fits.open(filename)
    header = sp[0].header

    wcs = WCS(header)
    #make index array
    index = np.arange(header['NAXIS1'])

    wavelength = wcs.wcs_pix2world(index[:,np.newaxis], 0)
    wavelength = wavelength.flatten()
    flux = sp[0].data

    date_obs = header['Date-OBS']
    return wavelength, flux, date_obs

    # This function uses the Doppler equivalency between wavelength and velocity
    
def wave2doppler(w, w0):
    import astropy.units as u
    w0_equiv = u.doppler_optical(w0)
    w_equiv = w.to(u.km/u.s, equivalencies=w0_equiv)
    return w_equiv


#if __name__ == "__main__":
#    import sys
#    read_spec(sys.argv[1])


#if __name__ == "__main__":
#    import sys
#    region_around_line(sys.argv[1],sys.argv[2],sys.argv[3])


