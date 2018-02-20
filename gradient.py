from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

"""
Determines the gradient of B1E5 from the ammonia data. Prints the gradient in the x-direction and y-direction, the magnitude of the gradient, and its angle with respect to north, east of north being a positive angle and west of north being negative.
"""

#Opens fits image, grabs data and header, then closes.
f = fits.open("B1E5_NH3_final2.fits")
data = f[0].data
header = f[0].header
f.close()
# data.shape = (1,1814,15,14)
data = data[::,::,3:12,3:12]

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
        result += intensity[i] * np.exp(-4*np.log(2)*((v-v0-velocity[i])/dv)**2)
    return result

#Fits the ammonia function to spectral data.
def fit(spectra, vel):
    low = (0,0,0,0)
    high = (20,50,10,5)
    popt, pcov = curve_fit(nh3_spectrum, vel, spectra, p0=[1.8,1.5,7,0.5], bounds=(low,high))
    return nh3_spectrum(vel, *popt)

#Input: column index. Output: weighted average spectra over the entire column at that row.
def avg_row(col,vel):
    a, b, c = data[0].shape
    result = np.zeros(a)
    weight_sum = 0
    for i in range(b):
        spectra = data[0,::-1,i,col]
        if not np.isnan(spectra).any():
            weight = 1/np.std(spectra[0:20])
            result += spectra*weight
            weight_sum += weight
    if weight_sum==0:
        return np.array([np.nan for i in range(a)])
    result /= weight_sum
    return fit(result,vel)

#Input: row index. Output: weighted average spectra over the entire column at that row.
def avg_col(row,vel):
    a, b, c = data[0].shape
    result = np.zeros(a)
    weight_sum = 0
    for i in range(c):
        spectra = data[0,::-1,row,i]
        if not np.isnan(spectra).any():
            weight = 1/np.std(spectra[0:20])
            result += spectra*weight
            weight_sum += weight
    if weight_sum==0:
        return np.array([np.nan for i in range(a)])
    result /= weight_sum
    return fit(result,vel)

#Input: list of average spectra and velocities. Output: list of average velocities
def get_avg_vels(vel, spec_list):
    result = []
    for spec in spec_list:
        weight_sum = 0
        avg = 0
        for i in range(len(spec)):
            if 6 < vel[i] < 8.5:
                weight = spec[i]
                avg += vel[i]*weight
                weight_sum += weight
        result.append(avg/weight_sum)
    print result
    return np.array(result)

#Returns a bunch of average rows that will be used to construct the intensities for the RA contour
def list_cols(vel):
    a, b, c = data[0].shape
    result = []
    for i in range(c):
        result.append(avg_row(i,vel))
    return np.array(result)

#Returns a bunch of average columns that will be used to construct the intensities for the dec contour
def list_rows(vel):
    a, b, c = data[0].shape
    result = []
    for i in range(b):
        result.append(avg_col(i,vel))
    return np.array(result)

#Gets values for the RA
def get_RA():
    a, b, c = data[0].shape
    pix = header['CRPIX1']
    delta = header['CDELT1']
    val = header['RA']
    RA = (np.arange(c)-pix)*delta + val
    return RA

#Gets values for the dec
def get_dec():
    a, b, c = data[0].shape
    pix = header['CRPIX2']
    delta = header['CDELT2']
    val = header['DEC']
    dec = (np.arange(b)-pix)*delta + val
    return dec

#Plots a contour of the right ascension and velocity, returning the slope of the best fit line (the gradient in the x-direction)
def RA_contour():
    res = header['CDELT3']/-1000
    v = 777
    a, b, c = data[0].shape
    first = -v*res
    last = (a-v)*res
    vel = np.arange(first,last,res)
    z = list_cols(vel)
    RA = get_RA()
    x, y = np.meshgrid(vel, RA)
    plt.contour(y,x,z)
    p = np.polyfit(RA, get_avg_vels(vel,z),1)
    plt.plot(RA, RA*p[0]+p[1],"k-")
    plt.xlim(plt.xlim()[::-1])
    plt.savefig("writeup_pics/RA_contour.png", bbox_inches="tight")
    plt.show()
    return p[0]

#Plots a contour of the declination and velocity, returning the slope of the best fit line (the gradient in the y-direction)
def dec_contour():
    res = header['CDELT3']/-1000
    v = 777
    a, b, c = data[0].shape
    first = -v*res
    last = (a-v)*res
    vel = np.arange(first,last,res)
    z = list_rows(vel)
    dec = get_dec()
    x, y = np.meshgrid(vel, dec)
    plt.contour(y,x,z)
    p = np.polyfit(dec, get_avg_vels(vel,z),1)
    plt.plot(dec, dec*p[0]+p[1],"k-")
    plt.savefig("writeup_pics/dec_contour.png", bbox_inches="tight")
    plt.show()
    return p[0]

#gets the magnitude and the angle of the gradient
def get_gradient():
    RA = RA_contour()
    dec = dec_contour()
    RA /= 235*np.pi/180
    dec /= 235*np.pi/180
    magnitude = np.sqrt(RA*RA + dec*dec)
    theta = np.arctan(RA/dec)
    print RA, dec
    print magnitude, theta*180/np.pi

get_gradient()
