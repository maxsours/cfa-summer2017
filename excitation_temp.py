import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit

data_file = fits.open('B1E2_NH3_final.fits')
data = data_file[0].data
header = data_file[0].header
data_file.close()

c = 300000 #speed of light in km/s
k = 1.38064852e-23 #Boltzmann constant
h = 6.62607004e-34 #Plank's constant
nu = 23694506000 #frequency for middle line

def nh3_spectrum(v, a, t0, v0, dv):
    return a * (1 - np.exp(-tv(v,t0,v0,dv)))

intensity = [0.44,0.22,0.50,0.06,0.28,0.10,0.28,1.40,0.06,0.90,0.10,0.06,0.11,0.06,0.50,0.28,0.44,0.22]

frequency = [1545.1,1534.1,617.7,582.8,571.8,19.9,16.8,10.5,5.9,-15.1,-24.5,-25.5,-36.5,-581.0,-590.3,-623.3,-1526.8,-1569]

ref_freq = 23694495.5
nu = ref_freq*1000
rest_freq = header['RESTFREQ']/1000
velocity = []
for f in frequency:
    velocity.append(-c*(f + ref_freq - rest_freq)/rest_freq)

def tv(v, t0, v0, dv):
    result = 0
    for i in range(0,18):
        result += intensity[i] * np.exp(-4*np.log(2)*((v-v0-velocity[i])/dv)**2)
    return t0 * result

def J(T):
    return (h*nu/k)*1/(np.exp(h*nu/(k*T))-1)

def solve(func, target, bounds, num):
    lower, upper = bounds
    for n in np.linspace(lower,upper,num):
        if func(n) >= target:
            return n
    return np.nan
"""
res = header['CDELT3']/-1000
v = 777
a, b, c = data[0].shape
spectra = data[0,::-1,7,7]
first = -v*res
last = (a-v)*res
vel = np.arange(first,last,res)
popt, pcov = curve_fit(nh3_spectrum, vel, spectra, p0=[1,1,7,1], bounds=(0,100))
x = np.linspace(2.73, 30, 10000)
y = J(x)
plt.plot(x,y)
plt.show()
print popt[0]
print popt[1]
print J(2.73)
print popt[0]+J(2.73)
print solve(J, popt[0]+J(2.73), (2.73,120), 100000)
"""
