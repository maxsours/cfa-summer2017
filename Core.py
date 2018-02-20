import numpy as np

"""
This class creates the code for Core objects, which are used to keep track of core properties. This is especially useful for creating scatter plots.
"""

class Core:
    
    def __init__(self, tau, Sv, n_nh3, Tex, Tk, intensity, mass, density, loc, types):
        self.tau = tau # optical depth
        self.Sv = Sv # velocity dispersion in km/s
        self.n_nh3 = n_nh3 # log column density
        self.Tex = Tex # excitation temperature in K
        self.Tk = Tk # kinetic temperature in K
        self.intensity = intensity # integrated intensity in K km/s
        self.mass = mass # in solar masses
        self.density = density # log of H2/cm3
        self.loc = loc # perseus or pipe
        self.types = types # starless or protostellar
        # ratio of non-thermal velocity dispersion to the speed of sound in the core
        self.Snt_to_Cs = self.get_Snt_to_Cs()
        self.vir_param = self.get_vir_param()

    def __str__(self):
        return "[{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}, {11}]".format(self.tau, self.Sv, self.n_nh3, self.Tex, self.Tk, self.intensity, self.mass, self.density, self.loc, self.types, self.Snt_to_Cs, self.vir_param)

    def get_vir_param(self):
        if self.mass is not None and self.density is not None and self.Sv is not None:
            radius = 1e-5*(0.75*self.mass*1.9891e33/(4.68e-24*np.pi*10**self.density))**(1/3.0) # in km
            constant = 7.5327226e-12 # (km/s)^2 * km / (G * Msun)
            return 5*constant*(self.Sv)**2*radius/self.mass
        return None

    def get_Snt_to_Cs(self):
        if self.Tk is not None and self.Sv is not None:
            Snt = np.sqrt((self.Sv)**2 - 1.38064852e-29*self.Tk/(17.03*1.6737e-27))
            Cs = np.sqrt(1.38064852e-29*self.Tk/(2.8*1.6737e-27))
            result = Snt/Cs
            if np.isfinite(result):
                return result
        return None

    # returns an instance variable corresponding to a specific keyword
    def get(self, kwd):
        data = [self.Tk, self.tau, self.Tex, self.Sv, self.intensity, self.n_nh3, self.mass, self.density, self.Snt_to_Cs, self.vir_param, self.types]
        names = ["Tk", "tau", "Tex", "Sv", "intensity", "n_nh3", "mass", "density", "Snt_to_Cs", "vir_param", "types"]
        dic = dict(zip(names, data))
        return dic[kwd]

