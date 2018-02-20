from Core import Core

"""
CoreList objects contain many cores and are used to get arrays of single values for plotting.
"""

class CoreList:

    def __init__(self):
        self.cores = []

    # adds a core to the list
    def add_core(self, core):
        self.cores.append(core)

    # Precodition: all arrays are the same length. The last two parameters are NOT arrays.
    def add_from_array(self, tau, Sv, n_nh3, Tex, Tk, intensity, mass, density, loc, types):
        for i in range(len(tau)):
            new_core = Core(tau[i], Sv[i], n_nh3[i], Tex[i], Tk[i], intensity[i], mass[i], density[i], loc, types[i])
            self.add_core(new_core)

    # Input: A keyword corresponding to an instance variable. Output: An array of that instance variable from the cores in this list
    def one_array(self, kwd):
        result = []
        for core in self.cores:
            result.append(core.get(kwd))
        return result
    
    # Input: Two keywords corresponding to two separate instance variables. Output: Two arrays of those instance variable from the cores in this list.
    def two_arrays(self, k1, k2):
        ar1 = self.one_array(k1)
        ar2 = self.one_array(k2)
        result1 = []
        result2 = []
        for i in range(len(ar1)):
            if ar1[i] is not None and ar2[i] is not None:
                result1.append(ar1[i])
                result2.append(ar2[i])
        return result1, result2
    
    # returns a core specified by a numerical keyword - which the cores should not have. Hopefully this method will never be used again.
    def get_core(self, num):
        for core in self.cores:
            if core.types==num:
                return core
        return

    # Adds the contents of another core list to this one.
    def add_cores(self, other_core_list):
        result = self
        for core in other_core_list.cores:
            result.add_core(core)
        return result

    # Returns a CoreList of all cores of a single type: perseus, pipe, starless, or protostellar. Note that starless and protostellar cores are also perseus cores.
    def get_type(self, ctype):
        if ctype=="perseus" or ctype=="pipe":
            result = CoreList()
            for core in self.cores:
                if core.loc==ctype:
                    result.add_core(core)
            return result
        if ctype=="starless" or ctype=="protostellar":
            result = CoreList()
            for core in self.cores:
                if core.types==ctype:
                    result.add_core(core)
            return result
        return


