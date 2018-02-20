import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit

"""
Fits ammonia spectra, ploting the data and the curve fit. Can be used on a variety of fits images.
"""

#Opens fits image, grabs data and header, then closes.
data_file = fits.open('B1E2_NH3_final2.fits')
data = data_file[0].data
header = data_file[0].header
data_file.close()

#data.shape = (1,1814,15,14)

c = 300000 #speed of light in km/s

#Function to create the ammonia spectrum, including the hyperfine lines
def nh3_spectrum(v, a, t0, v0, dv):
    return a * (1 - np.exp(-t0*tv(v,v0,dv)))

#Frequencies and relative intensities of the hyperfine lines.
intensity = [0.44,0.22,0.50,0.06,0.28,0.10,0.28,1.40,0.06,0.90,0.10,0.06,0.11,0.06,0.50,0.28,0.44,0.22]
frequency = [1545.1,1534.1,617.7,582.8,571.8,19.9,16.8,10.5,5.9,-15.1,-24.5,-25.5,-36.5,-581.0,-590.3,-623.3,-1526.8,-1569]

#Get the velocities of the hyperfine lines from their frequencies.
ref_freq = 23694495.5
rest_freq = header['RESTFREQ']/1000
velocity = []
for f in frequency:
    velocity.append(-c*(f + ref_freq - rest_freq)/rest_freq)

#A separate function used to compute the function for the ammonia spectrum.
def tv(v, v0, dv):
    result = 0
    for i in range(0,18):
        #print str(intensity[i]) + " " + str(velocity[i])
        result += intensity[i] * np.exp(-4*np.log(2)*((v-v0-velocity[i])/dv)**2)
    return result

#Plots the data, fitting the spectra and plotting both the spectra and the fit. It prints out the value determined by the fitting of each parameter of nh3_spectrum. It also prints out the error of the spectra in relation to the fit.
def plot(vel, spectra):
    low = (0,0,0,0)
    high = (20,50,10,5)
    plt.plot(vel,spectra,"r--")
    popt, pcov = curve_fit(nh3_spectrum, vel, spectra, p0=[1.8,1.5,7,0.5], sigma=error, bounds=(low,high))
    print popt
    print np.sqrt(np.diag(pcov))
    plt.plot(vel, nh3_spectrum(vel, *popt), "k-")
    plt.show()

#Determines the correct velocity values, and the noise of the spectra, then calls the plot method
res = header['CDELT3']/-1000
v = 777
a, b, c = data[0].shape
spectra = data[0,::-1,7,7] #the last two numbers can be modified to fit the spectra of any pixel in the fits image.
first = -v*res
last = (a-v)*res
vel = np.arange(first,last,res)
error = np.std(spectra[0:20])
plot(vel, spectra)


