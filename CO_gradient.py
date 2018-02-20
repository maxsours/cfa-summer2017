from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

"""
Determines the gradient of the entire B1E cloud, using carbon monoxide data. Prints out magnitude and angle with respect to north, east of north being positive and west of north being negative.
"""

#Opens the file, takes the data and the header, then closes the file.
f = fits.open("B1E_13CO_10_36_tan_cen.fits")
data = f[0].data
header = f[0].header
f.close()

#Input: column index. Output: weighted average spectra over the entire column at that row.
def avg_row(col):
    a, b, c = data.shape
    result = np.zeros(a)
    count = 0
    for i in range(b):
        spectra = data[::-1,i,col]
        if not np.isnan(spectra).any():
            result += spectra
            count += 1
    if count==0:
        return None
    result /= count
    return result

#Input: row index. Output: weighted average spectra over the entire column at that row.
def avg_col(row):
    a, b, c = data.shape
    result = np.zeros(a)
    count = 0
    for i in range(c):
        spectra = data[::-1,row,i]
        if not np.isnan(spectra).any():
            result += spectra
            count += 1
    if count==0:
        return None
    result /= count
    return result

#Input: list of average spectra and velocities. Output: list of average velocities
def get_avg_vels(vel, spec_list):
    result = []
    for spec in spec_list:
        weight_sum = 0
        avg = 0
        i = len(vel)-1
        while vel[i] > 5 and i >= 0:
            weight = spec[i]
            avg += vel[i]*weight
            weight_sum += weight
            i -= 1
        result.append(avg/weight_sum)
    return np.array(result)

#Returns a bunch of average rows that will be used to construct the intensities for the RA contour
def list_cols():
    a, b, c = data.shape
    result = []
    for i in range(c):
        result.append(avg_row(i))
    return np.array(result)

#Returns a bunch of average columns that will be used to construct the intensities for the dec contour
def list_rows():
    a, b, c = data.shape
    result = []
    for i in range(b):
        result.append(avg_col(i))
    return np.array(result)

#Gets values for the RA
def get_RA():
    a, b, c = data.shape
    pix = header['CRPIX1']
    delta = header['CDELT1']
    val = header['CRVAL1']
    RA = (np.arange(c)-pix)*delta + val
    return RA

#Gets values for the dec
def get_dec():
    a, b, c = data.shape
    pix = header['CRPIX2']
    delta = header['CDELT2']
    val = header['CRVAL2']
    dec = (np.arange(b)-pix)*delta + val
    return dec

#plots contour plot for the right ascension and velocity, returns the slope of the gradient in the x-direction
def RA_contour():
    res = header['CDELT3']/1000
    vpix = header['CRPIX3']
    a, b, c = data.shape
    vref = header['CRVAL3']/1000
    first = -vpix*res + vref
    last = (a-vpix)*res + vref
    vel = np.arange(first,last,res)
    vel = vel[::-1]
    z = list_cols()
    RA = get_RA()
    x, y = np.meshgrid(vel, RA)
    plt.contour(y,x,z,levels=[0.5,1,2,3,4,5])
    p = np.polyfit(RA, get_avg_vels(vel,z),1)
    plt.plot(RA, RA*p[0]+p[1],"k-")
    plt.xlim(plt.xlim()[::-1])
    plt.show()
    return p[0]

#plots contour plot for declination and velocity, returns the slope of the gradient in the y-direction.
def dec_contour():
    res = header['CDELT3']/1000
    vpix = header['CRPIX3']
    a, b, c = data.shape
    vref = header['CRVAL3']/1000
    first = -vpix*res + vref
    last = (a-vpix)*res + vref
    vel = np.arange(first,last,res)
    vel = vel[::-1]
    z = list_rows()
    dec = get_dec()
    x, y = np.meshgrid(vel, dec)
    plt.contour(y,x,z,levels=[0.5,1,2,3,4,5])
    p = np.polyfit(dec, get_avg_vels(vel,z),1)
    plt.plot(dec, dec*p[0]+p[1],"k-")
    plt.show()
    return p[0]

#gets the magnitude and the angle 
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
