import numpy as np
import matplotlib.pyplot as plt
from pipe_data import PipeData, PerseusData, Starless, Protostellar, PerseusNew
import matplotlib.axes as axes
import matplotlib.patches as mpatches
from CoreList import CoreList
from Core import Core
"""
This program contains a number of methods dedicated to taking the core data in pipe_data.py and creating charts and scatter plots from it.
"""

# get virial parameter
B1E2 = {"tau": 1.316, "Sv": 0.3493/2.3548, "n_nh3": 14.14, "Tex": 4.864, "Tk": 10.3, "intensity": 0.524, "mass(min)": 0.6, "mass(max)": 1.4, "density(min)": np.log10(55354), "density(max)": np.log10(129160)}
B1E5 = {"tau": 1.088, "Sv": 0.926/2.3548, "n_nh3(min)": 14.06, "n_nh3(max)":14.34, "Tex": 3.3369, "Tk(min)": 10, "Tk(max)": 19, "intensity": 0.179, "mass(min)": 0.5, "mass(max)": 2.0, "density(min)": np.log10(18839), "density(max)": np.log10(75358)}

titles = {"tau": "Optical depth", "Sv": "Velocity dispersion (km/s)", "n_nh3": "Column density (log(NH3))", "Tex": "Excitation temperature (K)", "Tk": "Kinetic temperature (K)", "intensity": "Integrated intensity (K km/s)", "mass": "Mass (solar)", "density": "Density (log(n/cm3))"}

k = 1.38064852e-29 # km2 kg s-2 K-1
nh3 = 17.03 # mean molecular weight of ammonia
h2 = 2.8 # mean molecular weight of the gas
mh = 1.6737e-27 # in kg

# master core list
cores = PipeData.core_list.add_cores(PerseusNew.core_list)

# computes the ratio of non-thermal velocity dispersion to the speed of sound in the nebula. This function has been integrated into the CoreList class, so it should not be needed anymore.
def Snt_to_Cs(Sv, Tk):
    Snt = np.sqrt(Sv*Sv - k*Tk/(nh3*mh))
    Cs = np.sqrt(k*Tk/(h2*mh))
    return Snt/Cs

#Jeans mass of B1E in solar masses assuming a temperature of 14K and a number density of 3000.
def jeans_mass():
    #5.347e21/1.989e30 = (km/s)^3 * pi^5/2 / (6 * G^3/2 * (kg/cm^3)^0.5) * (Msun/kg)
    print 5.347e21/1.989e30 * (np.sqrt(k*20/(h2*mh)))**3/(3000*h2*mh)**0.5 # = ~6.00 Msun

# Computes the virial parameter of a core
def vir_param(Sv, mass, density):
    radius = 1e-5*(0.75*mass*1.9891e33/(4.68e-24*np.pi*10**density))**(1/3.0) # in km
    constant = 7.5327226e-12 # (km/s)^2 * km / (G * Msun)
    return 5*constant*Sv*Sv*radius/mass

# plots a histogram of the ratio of non-thermal velocity dispersion to sound speed in the dense core. <1 = subsonic, >1 = supersonic.
def plot_supersonic():
    data = clear_nones(cores.one_array("Snt_to_Cs"))
    plt.hist(np.array(data))
    plt.title("Ratio of non-thermal velocity dispersion to speed of sound - Protostellar")
    # plots dashed lines where the values of B1E2 (red) and B1E5 (magenta) should be. B1E5 has a range of values, so there are two magenta lines at the maximum and minimum values.
    E2 = plt.axvline(Snt_to_Cs(B1E2["Sv"], B1E2["Tk"]), color="r", linestyle="dashed")
    E5 = plt.axvline(Snt_to_Cs(B1E5["Sv"], B1E5["Tk(min)"]), color="#808080")
    plt.axvline(Snt_to_Cs(B1E5["Sv"], B1E5["Tk(max)"]), color="#808080")
    plt.legend((E2, E5), ("B1E2", "B1E5"), scatterpoints=1, loc="upper right", ncol=2, fontsize=8)
    plt.savefig("writeup_pics/ratio_Snt_Cs.png", bbox_inches="tight")
    plt.show()

#Plots virial parameter on a histogram
def plot_vir_param():
    data = clear_nones(cores.one_array("vir_param"))
    plt.hist(np.array(data))
    plt.title("Virial Parameter")
    E2 = plt.axvline(vir_param(B1E2["Sv"], B1E2["mass(min)"], B1E2["density(min)"]), color="r", linestyle="dashed")
    plt.axvline(vir_param(B1E2["Sv"], B1E2["mass(max)"], B1E2["density(max)"]), color="r", linestyle="dashed")
    E5 = plt.axvline(vir_param(B1E5["Sv"], B1E5["mass(min)"], B1E5["density(min)"]), color="#808080")
    plt.axvline(vir_param(B1E5["Sv"], B1E5["mass(max)"], B1E5["density(max)"]), color="#808080")
    plt.legend((E2, E5), ("B1E2", "B1E5"), scatterpoints=1, loc="upper right", ncol=2, fontsize=8)
    plt.savefig("writeup_pics/vir_param_histogram.png", bbox_inches="tight")
    plt.show()
    

# generates the median value for each data point that is not split into starless and protostellar cores
def median(kwd):
    print "Pipe:"
    print "{0}: {1}".format(kwd, np.median(clear_nones(cores.get_type("pipe").one_array(kwd))))
    print "{0}: {1}".format(kwd, np.median(clear_nones(PipeData.dic[kwd])))
    print "\nPerseus:"
    print "{0}: {1}".format(kwd, np.median(clear_nones(cores.get_type("perseus").one_array(kwd))))
    print "{0}: {1}".format(kwd, np.median(clear_nones(PerseusData.dic[kwd])))

# Returns the input list with all of the None values removed
def clear_nones(lis):
    result = []
    for l in lis:
        if l is not None:
            result.append(l)
    return result

# plots a histogram of data that are not split up into starless and protostellar cores
def plot(kwd):
    pidata = clear_nones(cores.get_type("pipe").one_array(kwd))
    plt.subplot(2,1,1)
    plt.hist(pidata)
    plt.title(titles[kwd]+" - Pipe Data")
    c1 = "r" #red
    c2 = "#808080" #grey
    # plots dashed lines where the values of B1E2 (red) and B1E5 (orange) should be. If B1E5 a range of values, there will be two orange lines at the maximum and minimum values.
    E2 = plt.axvline(B1E2[kwd], color=c1, linestyle="dashed")
    try:
        E5 = plt.axvline(B1E5[kwd], color=c2)
    except KeyError:
        E5 = plt.axvline(B1E5[kwd+"(min)"], color=c2)
        plt.axvline(B1E5[kwd+"(max)"], color=c2)
    plt.legend((E2, E5), ("B1E2", "B1E5"), scatterpoints=1, loc="upper right", ncol=2, fontsize=8)
    pedata = clear_nones(cores.get_type("perseus").one_array(kwd))
    plt.subplot(2,1,2)
    plt.hist(pedata)
    plt.title(titles[kwd]+" - Perseus Data")
    E2 = plt.axvline(B1E2[kwd], color=c1, linestyle="dashed")
    try:
        E5 = plt.axvline(B1E5[kwd], color=c2)
    except KeyError:
        E5 = plt.axvline(B1E5[kwd+"(min)"], color=c2)
        plt.axvline(B1E5[kwd+"(max)"], color=c2)
    plt.legend((E2, E5), ("B1E2", "B1E5"), scatterpoints=1, loc="upper right", ncol=2, fontsize=8)
    if kwd=="n_nh3":
        plt.legend((E2, E5), ("B1E2", "B1E5"), scatterpoints=1, loc="upper left", ncol=2, fontsize=8)
    plt.savefig("writeup_pics/"+kwd+"_histogram.png", bbox_inches="tight")
    plt.show()

# generates the median value of data separated into starless and protostellar.
def median_n(kwd):
    print "Pipe:"
    print "{0}: {1}".format(kwd, np.median(clear_nones(cores.get_type("pipe").one_array(kwd))))
    print "\nPerseus (starless):"
    print "{0}: {1}".format(kwd, np.median(clear_nones(cores.get_type("starless").one_array(kwd))))
    print "\nPerseus (protostellar):"
    print "{0}: {1}".format(kwd, np.median(clear_nones(cores.get_type("protostellar").one_array(kwd))))

# plots data that is separaed into starless and protostellar (mass and density).
def plot_n(kwd):
    pipedata = clear_nones(cores.get_type("pipe").one_array(kwd))
    plt.hist(pipedata)
    plt.title(titles[kwd]+" - Pipe Data")
    c1 = "r"
    c2 = "#808080" #grey
    # plots dashed lines where the values of B1E2 (red) and B1E5 (orange) should be. If B1E5 or B1E2 a range of values, there will be two corresponding lines at the maximum and minimum values.
    try:
        E2 = plt.axvline(B1E2[kwd], color=c1, linestyle="dashed")
    except KeyError:
        E2 = plt.axvline(B1E2[kwd+"(min)"], color=c1, linestyle="dashed")
        plt.axvline(B1E2[kwd+"(max)"], color=c1, linestyle="dashed")
    try:
        E5 = plt.axvline(B1E5[kwd], color=c2)
    except KeyError:
        E5 = plt.axvline(B1E5[kwd+"(min)"], color=c2)
        plt.axvline(B1E5[kwd+"(max)"], color=c2)
    plt.legend((E2, E5), ("B1E2", "B1E5"), scatterpoints=1, loc="upper right", ncol=2, fontsize=8)
    plt.savefig("writeup_pics/"+kwd+"_histogram_pipe.png", bbox_inches="tight")
    plt.show()
    starlessdata = clear_nones(cores.get_type("starless").one_array(kwd))
    plt.subplot(2,1,1)
    plt.hist(starlessdata, bins=30)
    plt.title(titles[kwd]+" - Perseus Starless Data")
    try:
        plt.axvline(B1E2[kwd], color=c1, linestyle="dashed")
    except KeyError:
        plt.axvline(B1E2[kwd+"(min)"], color=c1, linestyle="dashed")
        plt.axvline(B1E2[kwd+"(max)"], color=c1, linestyle="dashed")
    try:
        plt.axvline(B1E5[kwd], color=c2)
    except KeyError:
        plt.axvline(B1E5[kwd+"(min)"], color=c2)
        plt.axvline(B1E5[kwd+"(max)"], color=c2)
    plt.legend((E2, E5), ("B1E2", "B1E5"), scatterpoints=1, loc="upper right", ncol=2, fontsize=8)
    plt.subplot(2,1,2)
    protodata = clear_nones(cores.get_type("protostellar").one_array(kwd))
    n, bins, patches = plt.hist(protodata,bins=30)
    plt.title(titles[kwd]+" - Perseus Protostellar Data")
    try:
        plt.axvline(B1E2[kwd], color=c1, linestyle="dashed")
    except KeyError:
        plt.axvline(B1E2[kwd+"(min)"], color=c1, linestyle="dashed")
        plt.axvline(B1E2[kwd+"(max)"], color=c1, linestyle="dashed")
    try:
        plt.axvline(B1E5[kwd], color=c2)
    except KeyError:
        plt.axvline(B1E5[kwd+"(min)"], color=c2)
        plt.axvline(B1E5[kwd+"(max)"], color=c2)
    plt.legend((E2, E5), ("B1E2", "B1E5"), scatterpoints=1, loc="upper right", ncol=2, fontsize=8)
    plt.savefig("writeup_pics/"+kwd+"_histogram_perseus.png", bbox_inches="tight")
    plt.show()

def plot_density():
    kwd = "density"
    scuba = pipe_exclusion().one_array("density")
    scuba = clear_nones(list(scuba))
    scuba = scuba[-17:]
    result = CoreList()
    pipe_cores = cores.get_type("pipe")
    for core in pipe_cores.cores:
        if not(core.tau is None and core.Sv is None and core.n_nh3 is None and core.Tex is None and core.Tk is None and core.intensity is None and core.mass is not None and core.density is not None):
            result.add_core(core)
            #print core.density
    pipedata = result.one_array("density")
    plt.hist(scuba)
    sc = mpatches.Patch(color="b")
    pi = mpatches.Patch(color="c")
    plt.hist(pipedata, color="c")
    plt.title(titles[kwd]+" - Pipe Data")
    c1 = "r"
    c2 = "#808080" #grey
    # plots dashed lines where the values of B1E2 (red) and B1E5 (orange) should be. If B1E5 or B1E2 a range of values, there will be two corresponding lines at the maximum and minimum values.
    E2 = plt.axvline(B1E2[kwd+"(min)"], color=c1, linestyle="dashed")
    plt.axvline(B1E2[kwd+"(max)"], color=c1, linestyle="dashed")
    E5 = plt.axvline(B1E5[kwd+"(min)"], color=c2)
    plt.axvline(B1E5[kwd+"(max)"], color=c2)
    plt.legend((sc, pi, E2, E5), ("SCUBA data", "Pipe cores", "B1E2", "B1E5"), scatterpoints=1, loc="upper right", ncol=2, fontsize=8)
    plt.savefig("writeup_pics/"+kwd+"_histogram_pipe.png", bbox_inches="tight")
    plt.show()
    starlessdata = Starless.dic_n[kwd]
    plt.subplot(2,1,1)
    plt.hist(starlessdata, bins=30)
    plt.title(titles[kwd]+" - Perseus Starless Data")
    plt.axvline(B1E2[kwd+"(min)"], color=c1, linestyle="dashed")
    plt.axvline(B1E2[kwd+"(max)"], color=c1, linestyle="dashed")
    plt.axvline(B1E5[kwd+"(min)"], color=c2)
    plt.axvline(B1E5[kwd+"(max)"], color=c2)
    plt.legend((E2, E5), ("B1E2", "B1E5"), scatterpoints=1, loc="upper right", ncol=2, fontsize=8)
    plt.subplot(2,1,2)
    protodata = Protostellar.dic_n[kwd]
    n, bins, patches = plt.hist(protodata,bins=30)
    plt.title(titles[kwd]+" - Perseus Protostellar Data")
    plt.axvline(B1E2[kwd+"(min)"], color=c1, linestyle="dashed")
    plt.axvline(B1E2[kwd+"(max)"], color=c1, linestyle="dashed")
    plt.axvline(B1E5[kwd+"(min)"], color=c2)
    plt.axvline(B1E5[kwd+"(max)"], color=c2)
    plt.legend((E2, E5), ("B1E2", "B1E5"), scatterpoints=1, loc="upper right", ncol=2, fontsize=8)
    plt.savefig("writeup_pics/"+kwd+"_histogram_perseus.png", bbox_inches="tight")
    plt.show()

# creates a scatter plot of two different variables. Also generates red and orange stars for B1E2 and B1E5, respectively.
def scatter(k1, k2):
    if k1 is not "density" and k2 is not "density":
        x1, y1 = cores.get_type("pipe").two_arrays(k1, k2)
        pi = plt.scatter(np.array(x1),np.array(y1),c="c") #cyan
    x2, y2 = cores.get_type("starless").two_arrays(k1, k2)
    x3, y3 = cores.get_type("protostellar").two_arrays(k1, k2)
    color1 = "r" #red
    color2 = "#808080" #grey
    s = plt.scatter(np.array(x2),np.array(y2),c="b") #blue
    pr = plt.scatter(np.array(x3),np.array(y3),c="#ff8c00") #orange
    if B1E2.keys().count(k1) > 0 and B1E2.keys().count(k2) > 0:
        E2 = plt.scatter(B1E2[k1], B1E2[k2], s=100, c=color1, marker="*")
    elif B1E2.keys().count(k1) > 0:
        E2 = plt.scatter(B1E2[k1], B1E2[k2+"(min)"], s=100, c=color1, marker="*")
        plt.scatter(B1E2[k1], B1E2[k2+"(max)"], s=100, c=color1, marker="*")
        plt.plot(np.array([B1E2[k1], B1E2[k1]]), np.array([B1E2[k2+"(min)"], B1E2[k2+"(max)"]]), color=color1)
    elif B1E2.keys().count(k2) > 0:
        E2 = plt.scatter(B1E2[k1+"(min)"], B1E2[k2], s=100, c=color1, marker="*")
        plt.scatter(B1E2[k1+"(max)"], B1E2[k2], s=100, c=color1, marker="*")
        plt.plot(np.array([B1E2[k1+"(min)"], B1E2[k1+"(max)"]]), np.array([B1E2[k2], B1E2[k2]]), color=color1)
    else:
        E2 = plt.scatter(B1E2[k1+"(min)"], B1E2[k2+"(min)"], s=100, c=color1, marker="*")
        plt.scatter(B1E2[k1+"(max)"], B1E2[k2+"(max)"], s=100, c=color1, marker="*")
        plt.plot(np.array([B1E2[k1+"(min)"], B1E2[k1+"(max)"]]), np.array([B1E2[k2+"(min)"], B1E2[k2+"(max)"]]), color=color1)
    # B1E5
    if B1E5.keys().count(k1) > 0 and B1E5.keys().count(k2) > 0:
        E5 = plt.scatter(B1E5[k1], B1E5[k2], s=100, c=color2, marker="*")
    elif B1E5.keys().count(k1) > 0:
        E5 = plt.scatter(B1E5[k1], B1E5[k2+"(min)"], s=100, c=color2, marker="*")
        plt.scatter(B1E5[k1], B1E5[k2+"(max)"], s=100, c=color2, marker="*")
        plt.plot(np.array([B1E5[k1], B1E5[k1]]), np.array([B1E5[k2+"(min)"], B1E5[k2+"(max)"]]), color=color2)
    elif B1E5.keys().count(k2) > 0:
        E5 = plt.scatter(B1E5[k1+"(min)"], B1E5[k2], s=100, c=color2, marker="*")
        plt.scatter(B1E5[k1+"(max)"], B1E5[k2], s=100, c=color2, marker="*")
        plt.plot(np.array([B1E5[k1+"(min)"], B1E5[k1+"(max)"]]), np.array([B1E5[k2], B1E5[k2]]), color=color2)
    else:
        E5 = plt.scatter(B1E5[k1+"(min)"], B1E5[k2+"(min)"], s=100, c=color2, marker="*")
        plt.scatter(B1E5[k1+"(max)"], B1E5[k2+"(max)"], s=100, c=color2, marker="*")
        plt.plot(np.array([B1E5[k1+"(min)"], B1E5[k1+"(min)"]]), np.array([B1E5[k2+"(min)"], B1E5[k2+"(max)"]]), color=color2)
        plt.plot(np.array([B1E5[k1+"(max)"], B1E5[k1+"(max)"]]), np.array([B1E5[k2+"(min)"], B1E5[k2+"(max)"]]), color=color2)
        plt.plot(np.array([B1E5[k1+"(min)"], B1E5[k1+"(max)"]]), np.array([B1E5[k2+"(min)"], B1E5[k2+"(min)"]]), color=color2)
        plt.plot(np.array([B1E5[k1+"(min)"], B1E5[k1+"(max)"]]), np.array([B1E5[k2+"(max)"], B1E5[k2+"(max)"]]), color=color2)
    plt.title(titles[k1]+ " vs. " + titles[k2])
    plt.xlabel(titles[k1])
    plt.ylabel(titles[k2])
    try:
        plt.legend((pi, s, pr, E2, E5), ("Pipe", "Perseus - Starless", "Perseus - Protostellar", "B1E2", "B1E5"), scatterpoints=1, loc="upper right", ncol=2, fontsize=8)
        if k2=="Tex":
            plt.legend((pi, s, pr, E2, E5), ("Pipe", "Perseus - Starless", "Perseus - Protostellar", "B1E2", "B1E5"), scatterpoints=1, loc="lower right", ncol=2, fontsize=8)
        if k1=="Tex":
            plt.legend((pi, s, pr, E2, E5), ("Pipe", "Perseus - Starless", "Perseus - Protostellar", "B1E2", "B1E5"), scatterpoints=1, loc="upper left", ncol=2, fontsize=8)
    except UnboundLocalError:
        plt.legend((s, pr, E2, E5), ("Perseus - Starless", "Perseus - Protostellar", "B1E2", "B1E5"), scatterpoints=1, loc="upper right", ncol=2, fontsize=8)
    plt.savefig("writeup_pics/"+k1+"_vs_"+k2+".png", bbox_inches="tight")
    plt.show()

def pipe_exclusion():
    result = CoreList()
    pipe_cores = cores.get_type("pipe")
    for core in pipe_cores.cores:
        if core.tau is None and core.Sv is None and core.n_nh3 is None and core.Tex is None and core.Tk is None and core.intensity is None and core.mass is not None and core.density is not None:
            result.add_core(core)
    return result

def scatter_mass_density():
    k1 = "mass"
    k2 = "density"
    x1 = PipeData.mass[-17:]
    y1 = PipeData.density[-17:]
    x2, y2 = cores.get_type("starless").two_arrays(k1, k2)
    x3, y3 = cores.get_type("protostellar").two_arrays(k1, k2)
    c1 = "r" #red
    c2 = "#808080" #grey
    pi = plt.scatter(np.array(x1),np.array(y1),c="c")
    s = plt.scatter(np.array(x2),np.array(y2),c="b")
    pr = plt.scatter(np.array(x3),np.array(y3),c="#ff8c00")
    E2 = plt.scatter(B1E2[k1+"(min)"], B1E2[k2+"(min)"], s=100, c=c1, marker="*")
    plt.scatter(B1E2[k1+"(max)"], B1E2[k2+"(max)"], s=100, c=c1, marker="*")
    plt.plot(np.array([B1E2[k1+"(min)"], B1E2[k1+"(max)"]]), np.array([B1E2[k2+"(min)"], B1E2[k2+"(max)"]]), color=c1)
    E5 = plt.scatter(B1E5[k1+"(min)"], B1E5[k2+"(min)"], s=100, c=c2, marker="*")
    plt.scatter(B1E5[k1+"(max)"], B1E5[k2+"(max)"], s=100, c=c2, marker="*")
    plt.plot(np.array([B1E5[k1+"(min)"], B1E5[k1+"(max)"]]), np.array([B1E5[k2+"(min)"], B1E5[k2+"(max)"]]), color=c2)
    plt.title(titles[k1]+ " vs. " + titles[k2])
    plt.xlabel(titles[k1])
    plt.ylabel(titles[k2])
    plt.legend((pi, s, pr, E2, E5), ("Pipe", "Perseus - Starless", "Perseus - Protostellar", "B1E2", "B1E5"), loc="upper right", ncol=2, fontsize=8)
    plt.savefig("writeup_pics/"+k1+"_vs_"+k2+".png", bbox_inches="tight")
    plt.show()

def scatter_all():
    keywords = ["tau", "Sv", "n_nh3", "Tex", "Tk", "intensity", "mass", "density"]
    for i in range(len(keywords)):
        for j in range(i+1, len(keywords)):
            k1 = keywords[i]
            k2 = keywords[j]
            if not (k1 is "mass" and k2 is "density"):
                scatter(k1, k2)
            else:
                scatter_mass_density()

def plot_all():
    keywords = ["tau", "Sv", "n_nh3", "Tex", "Tk", "intensity", "mass", "density"]
    for kwd in keywords:
        if kwd=="mass":
            plot_n(kwd)
        elif kwd=="density":
            plot_density()
        else:
            plot(kwd)
    plot_supersonic()
    plot_vir_param()

#Input: number density. Output: free-fall time in millions of years.
def Tff(density):
    density = 4.68e-24*10**density
    return np.sqrt(3*np.pi/(32*6.67408e-11*density))/31556926/1e6

#Input: velocity dispersion in km/s, number density. Output: Mass needed to make the virial parameter equal to 1.
def solve_mass(Sv, density):
    for n in np.linspace(0.0001,1000,10000):
        if vir_param(Sv, n, density) <= 1:
            return n
    return np.nan

#Input: mass (solar masses), number density. Output: accretion rate in solar masses per million years. Precondition: the virial parameter of the object >= 1.
def acc(Sv, mass, density):
    new_mass = solve_mass(Sv, density)
    dmass = new_mass - mass
    return dmass/Tff(density)

#Input: mass in solar masses, number density. Output: escape velocity in km/s.
def Vesc(mass, density):
    radius = 1e-2*(0.75*mass*1.9891e33/(4.68e-24*np.pi*10**density))**(1/3.0) # in m
    return 1e-3*np.sqrt(2*mass*1.9891e30*6.67408e-11/radius)

#Input: velocity dispersion in km/s, mass in solar masses, number density. Output: ratio of velocity dispersion to escape velocity.
def Sv_to_Vesc(Sv, mass, density):
    return Sv/Vesc(mass, density)

def plot_Tff():
    densities = clear_nones(cores.one_array("density"))
    data = []
    for v in densities:
        data.append(Tff(v))
    plt.title("Free-fall time in millions of years")
    plt.hist(np.array(data))
    E2 = plt.axvline(Tff(B1E2["density(min)"]), color="r", linestyle="dashed")
    plt.axvline(Tff(B1E2["density(max)"]), color="r", linestyle="dashed")
    E5 = plt.axvline(Tff(B1E5["density(min)"]), color="#808080")
    plt.axvline(Tff(B1E5["density(max)"]), color="#808080")
    plt.legend((E2, E5), ("B1E2", "B1E5"), scatterpoints=1, loc="upper right", ncol=2, fontsize=8)
    plt.savefig("writeup_pics/Tff_histogram.png", bbox_inches="tight")
    plt.show()

def plot_acc():
    vel_disps = cores.one_array("Sv")
    masses = cores.one_array("mass")
    densities = cores.one_array("density")
    data = []
    for i in range(len(densities)):
        try:
            if acc(vel_disps[i], masses[i], densities[i]) > 0:
                data.append(acc(vel_disps[i], masses[i], densities[i]))
        except TypeError:
            pass
    plt.title("Accretion needed to get virial paramater <1 in solar masses per million years")
    plt.hist(np.array(data), bins=40)
    E2 = plt.axvline(acc(B1E2["Sv"], B1E2["mass(min)"], B1E2["density(min)"]), color="r", linestyle="dashed")
    E5 = plt.axvline(acc(B1E5["Sv"], B1E5["mass(min)"], B1E5["density(min)"]), color="#808080")
    plt.axvline(acc(B1E5["Sv"], B1E5["mass(max)"], B1E5["density(max)"]), color="#808080")
    plt.legend((E2, E5), ("B1E2", "B1E5"), scatterpoints=1, loc="upper right", ncol=2, fontsize=8)
    plt.savefig("writeup_pics/acc_histogram.png", bbox_inches="tight")
    plt.show()

def plot_Vesc():
    masses = cores.one_array("mass")
    densities = cores.one_array("density")
    data = []
    for i in range(len(densities)):
        try:
            if Vesc(masses[i], densities[i]) > 0:
                data.append(Vesc(masses[i], densities[i]))
        except TypeError:
            pass
    plt.title("Escape velocity in km/s")
    plt.hist(np.array(data))
    E2 = plt.axvline(Vesc(B1E2["mass(min)"], B1E2["density(min)"]), color="r", linestyle="dashed")
    plt.axvline(Vesc(B1E2["mass(max)"], B1E2["density(max)"]), color="r", linestyle="dashed")
    E5 = plt.axvline(Vesc(B1E5["mass(min)"], B1E5["density(min)"]), color="#808080")
    plt.axvline(Vesc(B1E5["mass(max)"], B1E5["density(max)"]), color="#808080")
    plt.legend((E2, E5), ("B1E2", "B1E5"), scatterpoints=1, loc="upper right", ncol=2, fontsize=8)
    plt.savefig("writeup_pics/Vesc_histogram.png", bbox_inches="tight")
    plt.show()

def plot_Sv_to_Vesc():
    vel_disps = cores.one_array("Sv")
    masses = cores.one_array("mass")
    densities = cores.one_array("density")
    data = []
    for i in range(len(densities)):
        try:
            if Sv_to_Vesc(vel_disps[i], masses[i], densities[i]) > 0:
                data.append(Sv_to_Vesc(vel_disps[i], masses[i], densities[i]))
        except TypeError:
            pass
    plt.title("Ratio of velocity dispersion to escape velocity - Starless")
    plt.hist(np.array(data))
    E2 = plt.axvline(Sv_to_Vesc(B1E2["Sv"], B1E2["mass(min)"], B1E2["density(min)"]), color="r", linestyle="dashed")
    plt.axvline(Sv_to_Vesc(B1E2["Sv"], B1E2["mass(max)"], B1E2["density(max)"]), color="r", linestyle="dashed")
    E5 = plt.axvline(Sv_to_Vesc(B1E5["Sv"], B1E5["mass(min)"], B1E5["density(min)"]), color="#808080")
    plt.axvline(Sv_to_Vesc(B1E5["Sv"], B1E5["mass(max)"], B1E5["density(max)"]), color="#808080")
    plt.legend((E2, E5), ("B1E2", "B1E5"), scatterpoints=1, loc="upper right", ncol=2, fontsize=8)
    plt.savefig("writeup_pics/Sv_to_Vesc_histogram.png", bbox_inches="tight")
    plt.show()

plot_Sv_to_Vesc()
