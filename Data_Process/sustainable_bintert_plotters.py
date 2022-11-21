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
        filenameH = glob.glob('data/Hermite/particle_trajectory/*')
        filenameGRX = glob.glob('data/GRX/particle_trajectory/*')
        filename = [natsort.natsorted(filenameH), natsort.natsorted(filenameGRX)]

        self.pop = [[ ], [ ]]
        self.dedt = [[ ], [ ]]
        self.dadt = [[ ], [ ]]
        self.bin_life = [[ ], [ ]]
        self.sys_bin = [[ ], [ ]]
        self.ter_life = [[ ], [ ]]
        self.sys_ter = [[ ], [ ]]
        
        for int_ in range(2):
            for file_ in range(len(filename[int_])):
                with open(filename[int_][file_], 'rb') as input_file:
                    print('Reading file', file_, ':', input_file)
                    data = pkl.load(input_file)

                    self.pop[int_].append(np.shape(data)[0])
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

                    self.sys_bin[int_].append(np.sum(bin_formed))
                    self.sys_ter[int_].append(np.sum(ter_formed))

                    self.bin_life[int_].append(np.mean(bin_life))
                    self.ter_life[int_].append(np.mean(ter_life))

                    self.dedt[int_].append(np.mean(temp_dedt))
                    self.dadt[int_].append(np.mean(temp_dadt))

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

        norm_min = (min(np.nanmin(self.bin_life[0]), np.nanmin(self.ter_life[0]),
                        np.nanmin(self.bin_life[1]), np.nanmin(self.ter_life[1])))
        norm_max = (max(np.nanmax(self.bin_life[0]), np.nanmax(self.ter_life[0]),
                        np.nanmax(self.bin_life[1]), np.nanmax(self.ter_life[1])))
        normalise_p1 = plt.Normalize((norm_min), (norm_max))
        p1ymax = np.log10(norm_max)

        p2ymin = np.log10(10**6*min(np.nanmin(self.dedt[0]), np.nanmin(self.dedt[1])))
        p2ymax = np.log10(10**6*max(np.nanmax(self.dedt[0]), np.nanmax(self.dedt[1])))
        p2xmin = np.log10(10**6*min(np.nanmin(self.dadt[0]), np.nanmin(self.dadt[1])))
        p2xmax = np.log10(10**6*max(np.nanmax(self.dadt[0]), np.nanmax(self.dadt[1])))

        integrator = ['Hermite', 'GRX']

        for int_ in range(2):
            ini_pop = np.unique(self.pop[int_])
            bin_formed = np.empty(len(ini_pop))
            ter_formed = np.empty(len(ini_pop))
            bin_life = np.empty(len(ini_pop))
            ter_life = np.empty(len(ini_pop))

            iter = -1
            for pop_ in ini_pop:
                iter += 1
                idx = np.where(self.pop[int_] == pop_)[0]
                bin_formed[iter] = np.mean(self.sys_bin[int_][idx])
                ter_formed[iter] = np.mean(self.sys_ter[int_][idx])
                bin_life[iter] = np.mean(self.bin_life[int_][idx])
                ter_life[iter] = np.mean(self.ter_life[int_][idx])
            bin_formed = np.asarray(bin_formed)
            ter_formed = np.asarray(ter_formed)
            bin_life = np.asarray(bin_life)
            ter_life = np.asarray(ter_life)

            fig, ax = plt.subplots()
            plot_ini.tickers_pop(ax, ini_pop)
            ax.set_title(integrator[int_])
            ax.set_xlabel(r'IMBH Population [$N$]')
            ax.set_ylabel(r'$\log_{10}\langle N_{sys} \rangle$')
            plt.scatter(ini_pop[bin_formed>0], np.log10(bin_formed[bin_formed>0]), edgecolors  = 'black', c = (bin_life[bin_formed>0]), norm = (normalise_p1), label = 'Stable Binary')
            plt.scatter(ini_pop[ter_formed>0], np.log10(ter_formed[ter_formed>0]), edgecolors  = 'black', c = (ter_life[ter_formed>0]), norm = (normalise_p1), marker = 's', label = 'Stable Triple')
            plt.colorbar().set_label(r'$\langle t_{sys} \rangle$ [kyr]')
            ax.set_ylim(0, 1.1*p1ymax)
            ax.legend()
            plt.savefig('figures/binary_hierarchical/sys_formation_N_plot'+str(integrator[int_])+'.pdf', dpi=300, bbox_inches='tight')

            fig = plt.figure(figsize=(12.5, 6))
            ax1 = fig.add_subplot(121)
            ax2 = fig.add_subplot(122)
            ax.set_title(integrator[int_])
            for ax_ in [ax1, ax2]:
                plot_ini.tickers(ax_, 'plot')
                ax_.yaxis.set_major_formatter(mtick_formatter)            
            ax1.set_xlim(0.9*p2xmin, 1.1*p2xmax)
            ax1.set_ylim(0.9*p2ymin, 1.1*p2ymax)        
            ax2.set_xlim(-1.1*10**(p2xmin), 1.1*10**(p2xmax))
            ax2.set_ylim(0.9*10**(p2ymin), 1.1*10**(p2ymax))
            ax1.set_ylabel(r'$\log_{10}\langle \dot{(1-e)} \rangle_{SMBH}$ [Gyr$^{-1}$]')
            ax2.set_ylabel(r'$\langle \dot{(1-e)} \rangle_{SMBH}$ [Gyr$^{-1}$]')
            ax1.set_xlabel(r'$\log_{10}\langle \dot{a}^{-1} \rangle_{SMBH}$ [pc$^{-1}$Gyr$^{-1}$]')
            ax2.set_xlabel(r'$\langle \dot{a}^{-1} \rangle_{SMBH}$ [pc$^{-1}$Gyr$^{-1}$]')
            ax1.xaxis.set_major_formatter(mtick_formatter)
            colour_axes = ax1.scatter(np.log10(10**6*self.dadt[int_]), np.log10(10**6*self.dedt[int_]), edgecolors = 'black', c = self.pop[int_])
            colour_axes = ax2.scatter((10**6*self.dadt[int_]), (10**6*self.dedt[int_]), edgecolors = 'black', c = self.pop[int_])

            plt.colorbar(colour_axes, ax = ax2, label = r'Initial Population')
            plt.savefig('figures/binary_hierarchical/dadt_dedt_plot'+str(integrator[int_])+'.pdf', dpi=300, bbox_inches='tight')

######### PLOT ALL NN yr vs. DISTANCE

cst = sustainable_sys()
cst.system_formation()