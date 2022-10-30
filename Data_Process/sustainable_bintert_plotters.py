from time import time
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
        self.stab_bin = []
        self.stab_ter = []
        self.dedt = []
        self.dadt = []
        
        for file_ in range(len(filename)):
            with open(filename[file_], 'rb') as input_file:
                data = pkl.load(input_file)
                self.pop.append(np.shape(data)[0])
                bin_val = 0 
                ter_val = 0
                temp_dedt = []
                temp_dadt = []

                for parti_ in range(np.shape(data)[0]): #Iterate over all particles present
                    if parti_ == 0:
                        pass
                    else:
                        temp_dedt.append((data.iloc[parti_][-2][8][0] - data.iloc[parti_][2][8][0])/(np.shape(data)[1]-1))
                        temp_dadt.append((data.iloc[parti_][-2][7][0] - data.iloc[parti_][2][7][0]).value_in(units.pc)/(np.shape(data)[1]-1))
                        if isinstance(data.iloc[parti_][-1][0], np.uint64):
                            for col_ in range(np.shape(data)[1]): #Iterate over the data time-steps
                                if data.iloc[parti_][col_][6][1] == data.iloc[parti_][col_-1][6][1] \
                                    and data.iloc[parti_][col_][6][1] == data.iloc[parti_][col_-2][6][1] \
                                        and abs(data.iloc[parti_][col_][8][1]) < 1 \
                                            and data.iloc[parti_][col_][7][1] < 0.15 | units.parsec:
                                    bin_val += 1

                                    if data.iloc[parti_][col_][6][2] == data.iloc[parti_][col_-1][6][2] \
                                        and data.iloc[parti_][col_][6][2] == data.iloc[parti_][col_-2][6][2] \
                                            and abs(data.iloc[parti_][col_][8][2]) < 1 \
                                                and data.iloc[parti_][col_][7][2] < 0.5 | units.parsec:
                                        ter_val += 1
                        else:
                            for col_ in range(np.shape(data)[1]):
                                if col_ < 3 or col_ == np.shape(data)[1]-1:
                                    pass
                                else:
                                    if data.iloc[parti_][col_][6][1] == data.iloc[parti_][col_-1][6][1] \
                                        and data.iloc[parti_][col_][6][1] == data.iloc[parti_][col_-2][6][1] \
                                            and abs(data.iloc[parti_][col_][8][1]) < 1 \
                                                and data.iloc[parti_][col_][7][1] < 0.15 | units.parsec:
                                        bin_val += 1

                                        if data.iloc[parti_][col_][6][2] == data.iloc[parti_][col_-1][6][2] \
                                            and data.iloc[parti_][col_][6][2] == data.iloc[parti_][col_-2][6][2] \
                                                and abs(data.iloc[parti_][col_][8][2]) < 1 \
                                                    and data.iloc[parti_][col_][7][2] < 0.5 | units.parsec:
                                            ter_val += 1

                self.stab_bin.append(bin_val/2)
                self.stab_ter.append(ter_val/3)
                self.dedt.append(np.mean(temp_dedt))
                self.dadt.append(np.mean(temp_dadt))

        self.dadt /= max(self.dadt) 
        self.dedt /= max(self.dadt)

        self.stab_bin = np.asarray(self.stab_bin)
        self.stab_ter = np.asarray(self.stab_ter)
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
        xints = [i for i in range(1+int(max(self.pop))) if i %10 == 0]
        mtick_formatter = mtick.FormatStrFormatter('%0.2f')

        ini_pop = np.unique(self.pop)
        bin_formed = np.empty(len(ini_pop))
        ter_formed = np.empty(len(ini_pop))

        iter = -1
        for pop_ in ini_pop:
            iter += 1
            idx = np.where(self.pop == pop_)[0]
            bin_formed[iter] = np.mean(self.stab_bin[idx])
            ter_formed[iter] = np.mean(self.stab_ter[idx])
        bin_formed = np.asarray(bin_formed)
        ter_formed = np.asarray(ter_formed)

        fig, ax = plt.subplots()
        plot_ini.tickers(ax)
        ax.set_xlabel(r'IMBH Population [$N$]')
        ax.set_ylabel(r'$N_{sys}$]')
        ax.scatter(ini_pop, bin_formed, edgecolors  = 'black', color = 'red', label = 'Binary')
        ax.scatter(ini_pop, ter_formed, edgecolors  = 'black', color = 'red', marker = 's', label = 'Triple')
        ax.set_xlim(5,105)
        ax.set_xticks(xints)
        ax.legend()
        plt.savefig('figures/sys_formation_N_plot.pdf', dpi=300, bbox_inches='tight')

        fig, ax = plt.subplots()
        plot_ini.tickers(ax)
        ax.set_title('Hardening Rate vs. Eccentricity Evolution')
        ax.set_xlabel(r'$\langle \dot{a} \rangle / \langle \dot{a} \rangle_{max} $ [pc/kyr]')
        ax.set_ylabel(r'$\langle \dot{e} \rangle / \langle \dot{a} \rangle_{max} $ [kyr$^{-1}$]')
        colour_axes = ax.scatter(self.dadt, self.dedt, coloredges = 'black', c = self.pop)
        plt.colorbar(colour_axes, ax = ax, label = r'Initial Population')
        #ax.set_xlim(1.05*min(self.dadt), 1.05*max(self.dadt))
        #ax.set_ylim(0, 1.05*max(self.dedt))
        ax.yaxis.set_major_formatter(mtick_formatter)
        ax.xaxis.set_major_formatter(mtick_formatter)
        plt.savefig('figures/dadt_dedt_plot.pdf', dpi=300, bbox_inches='tight')

cst = sustainable_sys()
cst.system_formation()