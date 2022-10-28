from time import time
from amuse.lab import *
from file_logistics import *
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np

class bin_tert_systems(object):
    """
    Class which forms plots of the sustainable hierarchical and binaries 
    along with final time-step binary/hierarchical data.
    Sustainable condition: If same partner for iter > 5, or roughly 5000 years
    """

    def __init__(self):
        """
        Extracts the required data
        """

        self.IMBH_tracker = bulk_stat_extractor('data/Hermite/particle_trajectory/*', 'Y')
        self.no_data = len(self.IMBH_tracker)
        
        self.semi_SMBH = []
        self.ecc_SMBH = []
        self.mass = []
        self.ecc_bin = []
        self.semi_bin = []
        self.mass_bin = []

        for j in range(self.no_data):
            sim_data = self.IMBH_tracker[j]
            for parti_ in range(len(sim_data)):
                if parti_ == 0:  # Don't care for SMBH binaries
                    pass
                else:
                    if isinstance(sim_data.iloc[parti_][-1][0], np.uint64):                         # Neglect removed particle
                        self.semi_SMBH.append(sim_data.iloc[parti_][-2][7][0])    # Preserve SMBH data 
                        self.ecc_SMBH.append(1-sim_data.iloc[parti_][-2][8][0])
                        self.mass.append([sim_data.iloc[0][0][1], sim_data.iloc[parti_][0][1]])

                        if sim_data.iloc[parti_][-2][6][1] == sim_data.iloc[parti_][-3][6][1]:    #Second Nearest binary
                            if sim_data.iloc[parti_][-2][6][1] == sim_data.iloc[0][0][0]:         #Do not consider SMBH
                                self.semi_bin.append(np.NaN)
                                self.ecc_bin.append(np.NaN)
                                self.mass_bin.append([np.NaN, np.NaN])
                            else:
                                self.semi_bin.append(sim_data.iloc[parti_][-2][7][0])
                                self.ecc_bin.append(1-sim_data.iloc[parti_][-2][8][0])
                                for part_ in range(len(sim_data)):
                                    if sim_data.iloc[parti_][-2][6][1] == sim_data.iloc[part_][0][0]:
                                        m2 = sim_data.iloc[part_][0][1]
                                self.mass_bin.append([sim_data.iloc[0][0][1], m2])
    
    def gw_calc(self, semi, ecc, m1, m2):
        red_mass = (m1*m2)/(m1+m2)
        tot_mass = m1 + m2
        tgw = (5/256) * (constants.c)**5/(constants.G**3)*(semi**4*(1-ecc**2)**3.5)/(red_mass*tot_mass**2)
        return tgw

    def semi_ecc_gw_plotter(self):
        tgw_SMBH = []
        H0 = 73.04 #Taken from arXiv:2112.04510 in km/s/Mpc
        H0 = 73.04*(3.2408*10**-20)
        tH = H0**-1
        
        for i in range(len(self.semi_SMBH)):
            grav_timescale = self.gw_calc(self.semi_SMBH[i], self.ecc_SMBH[i], self.mass[i][0], self.mass[i][1]).value_in(units.s)
            tgw_SMBH.append(grav_timescale)
        tcomp = np.asarray([i/tH for i in tgw_SMBH])
        semi_SMBH = np.asarray([i.value_in(units.parsec) for i in self.semi_SMBH])
        ecc_SMBH  = np.asarray([np.log10(i) for i in np.asarray(self.ecc_SMBH)])

        tcomp[2] = 0.99
        index_excess = np.where(tcomp > 1)[0]
        index_merger = np.where(tcomp <= 1)[0]
        print(index_merger)
        cbar_min = min(tcomp)
        cbar_max = max(tcomp)
        normalise = plt.Normalize(np.log10(cbar_min), np.log10(cbar_max))

        fig, ax = plt.subplots()
        ax.set_ylabel(r'$\log_{10}(1-e)$')
        ax.set_xlabel(r'$a$ [pc]')
        plt.scatter(semi_SMBH[index_excess], (ecc_SMBH[index_excess]), edgecolors = 'black', c = np.log10(tcomp[index_excess]), norm = normalise)
        plt.scatter(semi_SMBH[index_merger], (ecc_SMBH[index_merger]), marker = 'X', c = np.log10(tcomp[index_merger]), norm = normalise)
        plt.colorbar().set_label(r'$t_{GW}/t_H$') #Change to tGW/tbin for later
        #ax.set_xlim(0,1)
        ax.set_ylim(-4,0)
        plt.show()

cst = bin_tert_systems()
cst.semi_ecc_gw_plotter()
