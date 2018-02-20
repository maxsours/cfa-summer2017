import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit
from excitation_temp import J, solve

data_file = fits.open('B1E5_average.fits')
data = data_file[0].data
header = data_file[0].header
data_file.close()

c = 300000 #speed of light in km/s
k = 1.38064852e-23 #Boltzmann constant
h = 6.62607004e-34 #Plank's constant

def nh3_spectrum(v, a, t0, v0, dv):
    return a * (1 - np.exp(-t0*tv(v,v0,dv)))

intensity = [0.44,0.22,0.50,0.06,0.28,0.10,0.28,1.40,0.06,0.90,0.10,0.06,0.11,0.06,0.50,0.28,0.44,0.22]

frequency = [1545.1,1534.1,617.7,582.8,571.8,19.9,16.8,10.5,5.9,-15.1,-24.5,-25.5,-36.5,-581.0,-590.3,-623.3,-1526.8,-1569]

ref_freq = 23694495.5
rest_freq = header['RESTFREQ']/1000
velocity = []
for f in frequency:
    velocity.append(-c*(f + ref_freq - rest_freq)/rest_freq)

ref_freq = 23694495.5
rest_freq = header['RESTFREQ']
cdelt = header['CDELT1']
Nv = header['NAXIS1']
pixV = header['CRPIX1']
nu = (np.arange(Nv)-(pixV-1))*cdelt+rest_freq
v = c*(rest_freq - nu)/rest_freq

#print velocity

def tv(v, v0, dv):
    result = 0
    for i in range(0,18):
        #print str(intensity[i]) + " " + str(velocity[i])
        result += intensity[i]/sum(intensity) * np.exp(-4*np.log(2)*((v-v0-velocity[i])/dv)**2)
    return result

def nh3_1_1(Tex, Sv, t0):
    return 1.85e13*(1+np.exp(-1.137/Tex))/(1-np.exp(-1.137/Tex))*Sv/2.3548*t0

def S(j):
    if (not j==0) and j%3==0:
        return 2
    return 1

def Z():
    result = 0
    Tk = 10
    Trot = Tk/(1+(Tk/41.5)*np.log(1+0.6*np.exp(-15.6/Tk)))
    B = 298117*1000000
    C = 186726*1000000
    for j in range(5):
        result += (2*j+1)*S(j)*np.exp(-h*((B*j*(j+1)+(C-B)*j*j))/(k*Trot))
    return result

def Zj(j):
    Tk = 10
    Trot = Tk/(1+(Tk/41.5)*np.log(1+0.6*np.exp(-15.6/Tk)))
    B = 298117*1000000
    C = 186726*1000000
    return (2*j+1)*S(j)*np.exp(-h*((B*j*(j+1)+(C-B)*j*j))/(k*Trot))

def plot(vel, spectra):
    low = (0,0,0,0)
    #low = 0
    high = (20,50,10,5)
    #high = np.inf
    plt.plot(vel,spectra,"r--")
    popt, pcov = curve_fit(nh3_spectrum, vel, spectra, p0=[1.8,1.5,7,0.5], sigma=error, bounds=(low,high))
    print popt
    print np.sqrt(np.diag(pcov))
    Tex = solve(J, popt[0]+J(2.73), (2.73,120), 100000)
    print Tex
    t = nh3_1_1(Tex, popt[3], popt[1])
    print t
    total = t*Z()/Zj(1)
    print total
    plt.plot(vel, nh3_spectrum(vel, *popt), "k-")
    plt.show()


res = header['CDELT3']/-1000
spectra = data[0,0,0,::-1]

v = v[::-1]
error = spectra*0 + np.std(spectra[0:20])
plot(v, spectra)


