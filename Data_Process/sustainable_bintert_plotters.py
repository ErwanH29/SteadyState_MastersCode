from amuse.lab import *
from spatial_plotters import *
from file_logistics import *
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np

class sustainable_sys(object):
    """
    Class which extracts ALL information of particle_trajectory files in an memory-optimised manner
    and appends the required data values into arrays.
    """

    def __init__(self):
        filename = glob.glob('data/Hermite/particle_trajectory/*')
        filename = natsort.natsorted(filename)

        self.pop = []
        self.dedt = []
        self.dadt = []
        self.bin_life = []
        self.sys_bin = []
        self.ter_life = []
        self.sys_ter = []
        
        for file_ in range(len(filename)):
            with open(filename[file_], 'rb') as input_file:
                print('Reading file', file_, ':', filename[file_])
                data = pkl.load(input_file)

                self.pop.append(np.shape(data)[0])
                bin_key = []
                ter_key = []
                bin_life = []
                bin_formed = []
                ter_life = []
                ter_formed = []
                temp_dedt = []
                temp_dadt = []

                avg_mass = 0
                for parti_ in range(np.shape(data)[0]-1):
                    avg_mass +=  data.iloc[parti_+1][0][1].value_in(units.kg)
                avg_mass /= (np.shape(data)[0])

                for parti_ in range(np.shape(data)[0]):
                    bin_val = 0
                    ter_val = 0
                    if parti_ == 0:
                        pass
                    else:
                        temp_dedt.append((data.iloc[parti_][-2][8][0] - data.iloc[parti_][2][8][0])/(np.shape(data)[1]-3))
                        temp_dadt.append((data.iloc[parti_][-2][7][0]**-1 - data.iloc[parti_][2][7][0]**-1).value_in((units.pc)**-1)/(np.shape(data)[1]-3))
                        for col_ in range(np.shape(data)[1]-1):
                            if abs(data.iloc[parti_][col_][8][1]) < 1 \
                                and data.iloc[parti_][col_][7][1] < 0.02 | units.parsec:
                                mass1 = data.iloc[parti_][0][1]
                                for part_ in range(np.shape(data)[0]):
                                    if data.iloc[part_][0][0] == data.iloc[parti_][col_][6][1]:
                                        mass2 = data.iloc[part_][0][1]
                                        bin_BE = ((constants.G*mass1*mass2)/(data.iloc[parti_][col_][7][1])).value_in(units.J) 
                                        if bin_BE > (150000)**2*(1+mass1/mass2):   #Hard binary conditions based on Quinlan 1996b
                                            bin_val += 1
                                            bin_key.append(data.iloc[parti_][col_][6][1])

                                            #Calculate tertiary. The stability equality is based on Mardling and Aarseth 2001
                                            if data.iloc[parti_][col_][8][2] < 1 \
                                                and data.iloc[parti_][col_][7][2] < 0.1 | units.parsec:
                                                semi_inner = data.iloc[parti_][col_][7][1]
                                                semi_outer = data.iloc[parti_][col_][7][2]
                                                for part_ in range(np.shape(data)[0]):
                                                    if data.iloc[part_][0][0] == data.iloc[parti_][col_][6][2]:
                                                        mass_outer = data.iloc[part_][0][1]
                                                ecc_outer = data.iloc[parti_][col_][8][2]
                                                semi_ratio = semi_outer / semi_inner
                                                equality = 2.8 * ((1+mass_outer/(mass1+mass2))*(1+ecc_outer)/(1-ecc_outer**2)**0.5)**0.4
                                                if semi_ratio > equality:
                                                    ter_val += 1
                                                    ter_key.append(data.iloc[parti_][col_][6][2])

                                else:
                                        pass

                    bin_formed.append(len(np.unique(bin_key)))
                    ter_formed.append(len(np.unique(ter_key)))

                    if len(bin_formed) > 0:
                        bin_life.append(len(bin_key)/(max(1,bin_formed[-1])))
                    if len(ter_formed) > 0:
                        ter_life.append(len(ter_key)/(max(1,ter_formed[-1])))

                self.sys_bin.append(np.sum(bin_formed))
                self.sys_ter.append(np.sum(ter_formed))

                self.bin_life.append(np.mean(bin_life))
                self.ter_life.append(np.mean(ter_life))

                self.dedt.append(np.mean(temp_dedt))
                self.dadt.append(np.mean(temp_dadt))

        self.bin_life = np.asarray(self.bin_life)
        self.ter_life = np.asarray(self.ter_life)
        self.sys_bin = np.asarray(self.sys_bin)
        self.sys_ter = np.asarray(self.sys_ter)
        self.dedt = np.asarray(self.dedt)
        self.dadt = np.asarray(self.dadt)
        self.pop = np.asarray(self.pop)
        self.pop[self.pop %10 != 0] -= 1
        self.pop = np.asarray(self.pop)

    def system_formation(self):
        """
        Function to plot various 'sustainable system' plots
        """

        plot_ini = plotter_setup()
        mtick_formatter = mtick.FormatStrFormatter('%0.2f')

        ini_pop = np.unique(self.pop)
        bin_formed = np.empty(len(ini_pop))
        ter_formed = np.empty(len(ini_pop))
        bin_life = np.empty(len(ini_pop))
        ter_life = np.empty(len(ini_pop))

        iter = -1
        for pop_ in ini_pop:
            iter += 1
            idx = np.where(self.pop == pop_)[0]
            bin_formed[iter] = np.mean(self.sys_bin[idx])
            ter_formed[iter] = np.mean(self.sys_ter[idx])
            bin_life[iter] = np.mean(self.bin_life[idx])
            ter_life[iter] = np.mean(self.ter_life[idx])
        bin_formed = np.asarray(bin_formed)
        ter_formed = np.asarray(ter_formed)
        bin_life = np.asarray(bin_life)
        ter_life = np.asarray(ter_life)

        norm_min = (min(np.min(bin_life[bin_formed>0]), np.min(ter_life[ter_formed>0])))
        norm_max = (max(np.max(bin_life[bin_formed>0]), np.max(ter_life[ter_formed>0])))

        fig, ax = plt.subplots()
        plot_ini.tickers_pop(ax, ini_pop)
        ax.set_xlabel(r'IMBH Population [$N$]')
        ax.set_ylabel(r'$\log_{10}\langle N_{sys} \rangle$')
        normalise = plt.Normalize((norm_min), (norm_max))
        plt.scatter(ini_pop[bin_formed>0], np.log10(bin_formed[bin_formed>0]), edgecolors  = 'black', c = (bin_life[bin_formed>0]), norm = (normalise), label = 'Stable Binary')
        plt.scatter(ini_pop[ter_formed>0], np.log10(ter_formed[ter_formed>0]), edgecolors  = 'black', c = (ter_life[ter_formed>0]), norm = (normalise), marker = 's', label = 'Stable Triple')
        plt.colorbar().set_label(r'$\langle t_{sys} \rangle$ [kyr]')
        ax.set_ylim(0,1.1*np.log10(max(np.max(bin_formed), np.max(ter_formed))))
        ax.legend()
        plt.savefig('figures/binary_hierarchical/sys_formation_N_plot.pdf', dpi=300, bbox_inches='tight')

        fig = plt.figure(figsize=(12.5, 6))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        for ax_ in [ax1, ax2]:
            plot_ini.tickers(ax_, 'plot')
            ax_.yaxis.set_major_formatter(mtick_formatter)
        ax1.set_ylabel(r'$\log_{10}\langle \dot{(1-e)} \rangle_{SMBH}$ [Gyr$^{-1}$]')
        ax2.set_ylabel(r'$\langle \dot{(1-e)} \rangle_{SMBH}$ [Gyr$^{-1}$]')
        ax1.set_xlabel(r'$\log_{10}\langle \dot{a}^{-1} \rangle_{SMBH}$ [pc$^{-1}$Gyr$^{-1}$]')
        ax2.set_xlabel(r'$\langle \dot{a}^{-1} \rangle_{SMBH}$ [pc$^{-1}$Gyr$^{-1}$]')
        ax1.xaxis.set_major_formatter(mtick_formatter)
        colour_axes = ax1.scatter(np.log10(10**6*self.dadt), np.log10(10**6*self.dedt), edgecolors = 'black', c = self.pop)
        colour_axes = ax2.scatter((10**6*self.dadt), (10**6*self.dedt), edgecolors = 'black', c = self.pop)
        plt.colorbar(colour_axes, ax = ax2, label = r'Initial Population')
        plt.savefig('figures/binary_hierarchical/dadt_dedt_plot.pdf', dpi=300, bbox_inches='tight')

cst = sustainable_sys()
cst.system_formation()