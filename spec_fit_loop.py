import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit
from spectral_cube import SpectralCube

data_file = fits.open('B1E5_NH3_final_new.fits')
data = data_file[0].data
header = data_file[0].header
data_file.close()

c = 300000 #speed of light in km/s
k = 1.38064852e-23 #Boltzmann constant
T = 10.26 #for B1E2, temperature in K
nh3 = 17.03 #molecular weight of NH3
mh = 1.6737e-27 #mass of atomic hydrogen in kg

def nh3_spectrum(v, a, t0, v0, dv):
    return a * (1 - np.exp(-tv(v,t0,v0,dv)))

intensity = [0.44,0.22,0.50,0.06,0.28,0.10,0.28,1.40,0.06,0.90,0.10,0.06,0.11,0.06,0.50,0.28,0.44,0.22]

frequency = [1545.1,1534.1,617.7,582.8,571.8,19.9,16.8,10.5,5.9,-15.1,-24.5,-25.5,-36.5,-581.0,-590.3,-623.3,-1526.8,-1569]

ref_freq = 23694495.5
rest_freq = header['RESTFREQ']/1000
velocity = []
for f in frequency:
    velocity.append(-c*(f + ref_freq - rest_freq)/rest_freq)

res = header['CDELT3']/-1000
v = 777
a, b, c = data[0].shape
first = -v*res
last = (a-v)*res
vel = np.arange(first,last,res)

def tv(v, t0, v0, dv):
    result = 0
    for i in range(0,18):
        result += intensity[i]/sum(intensity) * np.exp(-4*np.log(2)*((v-v0-velocity[i])/dv)**2) # base e
    return t0 * result

f = fits.open("test4.fits", mode="update")
low = (0,0,0,0)
high = (20,50,10,5)
for i in range(b): # range(4,12)
    for j in range(c): # range(4,11)
        spectra = data[0,::-1,i,j]
        if not np.isnan(spectra).all():
            #"""            
            try:
                popt, pcov = curve_fit(nh3_spectrum, vel, spectra, p0=[1,1,7,1], bounds=(low,high))
                vel_disp = popt[3]/2.354820045
                f[0].data[i,j] = vel_disp
                #f[0].data[i,j] = popt[2]
            except RuntimeError:
                f[0].data[i,j] = np.nan
            """
            try:
                popt, pcov = curve_fit(nh3_spectrum, vel, spectra, p0=[1,1,7,1], bounds=(low,high))
                perr = np.sqrt(np.diag(pcov))
                if perr[2] < 0.1 and perr[3] < 0.3 and popt[3] > 0.07:
                    vel_disp = popt[3]/2.354820045 # velocity dispersion
                    St2=k*T/(1000000*nh3*mh) #square of thermal velocity dispersion
                    try:                
                        f[0].data[i,j] = np.sqrt(vel_disp**2 - St2)
                    except RuntimeWarning:
                        f[0].data[i,j] = np.nan
                    #f[0].data[i,j] = vel_disp
                else:
                    f[0].data[i,j] = np.nan
            except RuntimeError:
                f[0].data[i,j] = np.nan
                #"""
            
        else:
            f[0].data[i,j] = np.nan
f.close()

