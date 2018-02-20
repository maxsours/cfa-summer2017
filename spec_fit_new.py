# this code is dumb. Use spec_fit_average instead.

import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit

data_file = fits.open('B1E2_average.fits')
data = data_file[0].data
header = data_file[0].header
data_file.close()

c = 300000 #speed of light in km/s

def nh3_spectrum(v, a, t0, v0, dv):
    return a * (1 - np.exp(-t0*tv(v,v0,dv)))

intensity = [0.44,0.22,0.50,0.06,0.28,0.10,0.28,1.40,0.06,0.90,0.10,0.06,0.11,0.06,0.50,0.28,0.44,0.22]

frequency = [1545.1,1534.1,617.7,582.8,571.8,19.9,16.8,10.5,5.9,-15.1,-24.5,-25.5,-36.5,-581.0,-590.3,-623.3,-1526.8,-1569]

ref_freq = 23694495.5
rest_freq = header['RESTFREQ']/1000
velocity = []
for f in frequency:
    velocity.append(-c*(f + ref_freq - rest_freq)/rest_freq)

naxis1 = header['NAXIS1']
crpix1 = header['CRPIX1']
cdelt1 = header['CDELT1']
restfreq = header['RESTFREQ']

v_res = c*cdelt1/restfreq

nu = (np.arange(naxis1)-(crpix1-1))*cdelt1+restfreq

vel = c*(restfreq-nu)/restfreq

vel = vel[::-1]

def tv(v, v0, dv):
    result = 0
    for i in range(0,18):
        #print str(intensity[i]) + " " + str(velocity[i])
        result += intensity[i] * np.exp(-4*np.log(2)*((v-v0-velocity[i])/dv)**2)
    return result

def plot(vel, spectra):
    low = (0,0,0,0)
    #low = 0
    high = (20,50,10,5)
    #high = np.inf
    plt.plot(vel,spectra,"bo")
    popt, pcov = curve_fit(nh3_spectrum, vel, spectra, p0=[1.8,1.5,7,0.5], sigma=(error), bounds=(low,high))
    print popt
    print np.sqrt(np.diag(pcov))
    plt.plot(vel, nh3_spectrum(vel, *popt), "g^")
    #plt.show()

spectra = data[0,0,0,::-1]
error = spectra*0 + np.std(spectra[0:20])
plot(vel, spectra)


