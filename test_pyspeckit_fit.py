import matplotlib.pyplot as plt
plt.ion()

import pyspeckit

file_in='data/Core2_N2Dp_21_fine5.fits'
file_gauss= 'fits/Gaussian_fitted_parameters3.fits'
cube = pyspeckit.Cube(file_in)

F=False
T=True

cube.fiteach(fittype='gaussian',  guesses=[1.3, 6.7e3, 0.3e3], # # Peak=1.3, v_center=7.1, \sigma_v=0.3 km/s
             verbose_level=1, signal_cut=5, limitedmax=[F,T,F], maxpars=[0,8.0e3,1.0e3],
             limitedmin=[T,T,F], minpars=[0,5.0e3,0.05e3], use_neighbor_as_guess=True, 
             start_from_point=(10,10), multicore=4)

cube.write_fit(file_gauss, clobber=True)
cube.mapplot()
cube.plot_spectrum(10,10, plot_fit=True)
cube.mapplot.plane = cube.parcube[2,:,:]
cube.mapplot(estimator=None)


cube.load_model_fit(file_gauss, npars=3, npeaks=1)
# cube.mapplot()
# cube.plot_spectrum(10,10, plot_fit=True)
# cube.mapplot.plane = cube.parcube[2,:,:]
# cube.mapplot(estimator=None, vmin=0.05, vmax=0.3)
# plt.draw()
# plt.show()



import astropy.units as u
freq_line=154.217011*u.GHz
cube.xarr.refX = freq_line
cube.xarr.velocity_convention = 'radio'
cube.xarr.convert_to_unit('km/s')

from pyspeckit.spectrum.models import n2dp
cube.Registry.add_fitter('n2dp_vtau', pyspeckit.models.n2dp.n2dp_vtau_fitter, 4)

xmax=10; ymax=12
vmin=5.0; vmax=8.0
print('start optically thin fit')
cube.fiteach(fittype='n2dp_vtau',  guesses=[5, 0.1, 7.0, 0.1], # Tex=5K, tau=0.1, v_center=7.1, \sigma_v=0.3 km/s
             verbose_level=1, signal_cut=5,
             limitedmax=[F,F,T,T],
             limitedmin=[T,T,T,T],
             minpars=[ 0,  0,vmin,0.05],
             maxpars=[35.0,0,vmax,1.0],
             fixed=[F,T,F,F], 
             use_neighbor_as_guess=True, 
             start_from_point=(10,10), multicore=4)

cube.mapplot()
cube.plot_spectrum(xmax,ymax, plot_fit=True)
cube.mapplot.plane = cube.parcube[-1,:,:]
cube.mapplot(estimator=None, vmin=0.05, vmax=0.3)



# export files
parcube=cube.parcube
parcube[parcube==0]=np.nan
hd=cube.header
rm_key = ['NAXIS3', 'CRPIX3', 'CDELT3', 'CUNIT3', 'CTYPE3', 'CRVAL3', 'CROTA3']
rm_key = ['CROTA3']
for key_i in rm_key:
    hd.remove(key_i)
hd['NAXIS']=2
# Vlsr
hd['BUNIT']='km/s'
file_out='fits/2_thin_Vlsr.fits'
fits.writeto( file_out, parcube[2,::], hd, clobber=True)
# Sigma_v
hd['BUNIT']='km/s'
file_out='fits/2_thin_sigmav.fits'
fits.writeto( file_out, parcube[3,::], hd, clobber=True)