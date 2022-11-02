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
        filename = glob.glob('data/Hermite/particle_trajectory_temp/*')
        filename = natsort.natsorted(filename)

        self.pop = []
        self.stab_bin = []
        self.stab_ter = []
        self.stab_hard = []
        self.dedt = []
        self.dadt = []
        self.bin_life = []
        self.ter_life = []
        self.hard_life = []
        
        for file_ in range(len(filename)):
            with open(filename[file_], 'rb') as input_file:
                data = pkl.load(input_file)
                self.pop.append(np.shape(data)[0])
                bin_val = 0 
                ter_val = 0
                hard_val = 0
                bin_key = []
                ter_key = []
                hard_key = []
                temp_dedt = []
                temp_dadt = []

                avg_mass = 0
                for parti_ in range(np.shape(data)[0]-1):
                    avg_mass +=  data.iloc[parti_+1][0][1].value_in(units.kg)
                avg_mass /= (np.shape(data)[0])


                for parti_ in range(np.shape(data)[0]): #Iterate over all particles present
                    if parti_ == 0:
                        pass
                    else:
                        temp_dedt.append((data.iloc[parti_][-2][8][0] - data.iloc[parti_][2][8][0])/(np.shape(data)[1]-1))
                        temp_dadt.append((data.iloc[parti_][-2][7][0] - data.iloc[parti_][2][7][0]).value_in(units.pc)/(np.shape(data)[1]-1))
                        if isinstance(data.iloc[parti_][-1][0], np.uint64):
                            for col_ in range(np.shape(data)[1]): #Iterate over the data time-steps
                                if data.iloc[parti_][col_][6][1] == data.iloc[parti_][col_-1][6][1] \
                                    and abs(data.iloc[parti_][col_][8][1]) < 1 \
                                        and data.iloc[parti_][col_][7][1] < 0.15 | units.parsec:
                                    print('BIN')
                                    bin_val += 1
                                    bin_key.append(data.iloc[parti_][col_][6][1])

                                    #Calculate whether hard
                                    mass1 = data.iloc[parti_][0][1]
                                    for part_ in range(np.shape(data)[0]):
                                        if data.iloc[part_][0][0] == data.iloc[parti_][col_][6][1]:
                                            mass2 = data.iloc[part_][0][1]
                                            bin_BE = ((constants.G*mass1*mass2)/(2*data.iloc[parti_][col_][7][1])).value_in(units.J)
                                            if bin_BE > avg_mass * 150000:
                                                print('HArdBIN')
                                                hard_val += 1
                                                hard_key.append(data.iloc[parti_][col_][6][1])
                                        else:
                                            pass
                                    
                                    #Calculate tertiary
                                    if data.iloc[parti_][col_][6][2] == data.iloc[parti_][col_-1][6][2] \
                                        and abs(data.iloc[parti_][col_][8][2]) < 1 \
                                            and data.iloc[parti_][col_][7][2] < 0.5 | units.parsec:
                                        ter_val += 1
                                        ter_key.append(data.iloc[parti_][col_][6][2])
                        else:
                            for col_ in range(np.shape(data)[1]):
                                if col_ < 3 or col_ == np.shape(data)[1]-1:
                                    pass
                                else:
                                    if data.iloc[parti_][col_][6][1] == data.iloc[parti_][col_-1][6][1] \
                                        and abs(data.iloc[parti_][col_][8][1]) < 1 \
                                            and data.iloc[parti_][col_][7][1] < 0.15 | units.parsec:
                                        bin_val += 1
                                        bin_key.append(data.iloc[parti_][col_][6][1])
                                        print('BIN')

                                        #Calculate whether hard
                                        mass1 = data.iloc[parti_][0][1]
                                        for part_ in range(np.shape(data)[0]):
                                            if data.iloc[part_][0][0] == data.iloc[parti_][col_][6][1]:
                                                mass2 = data.iloc[part_][0][1]
                                                bin_BE = ((constants.G*mass1*mass2)/(2*data.iloc[parti_][col_][7][1])).value_in(units.J)
                                                if bin_BE > avg_mass * 150000:
                                                    hard_val += 1
                                                    hard_key.append(data.iloc[parti_][col_][6][1])
                                                    print('hard')
                                            else:
                                                pass

                                        if data.iloc[parti_][col_][6][2] == data.iloc[parti_][col_-1][6][2] \
                                                and abs(data.iloc[parti_][col_][8][2]) < 1 \
                                                    and data.iloc[parti_][col_][7][2] < 0.5 | units.parsec:
                                            ter_val += 1
                                            ter_key.append(data.iloc[parti_][col_][6][2])

                bin_val0 = bin_val
                hard_val0 = hard_val
                ter_val0 = ter_val
                for rep_ in range(len(bin_key)):
                    if bin_key[rep_] != 0:
                        bin_sust = -1
                        if rep_ == 0:
                            pass
                        else:
                            if bin_key[rep_] == bin_key[rep_-1]:
                                bin_val -= 1
                                bin_sust += 1

                for rep_ in range(len(ter_key)):
                    if ter_key[rep_] != 0:
                        ter_sust = -1
                        if rep_ == 0:
                            pass
                        else:
                            if ter_key[rep_] == ter_key[rep_-1]:
                                ter_val -= 1
                                ter_sust += 1
                
                for rep_ in range(len(hard_key)):
                    if hard_key[rep_] != 0:
                        hard_sust = -1
                        if rep_ == 0:
                            pass
                        else:
                            if hard_key[rep_] == hard_key[rep_-1]:
                                hard_val -= 1
                                hard_sust += 1
                
                self.stab_bin.append(bin_val)
                self.stab_ter.append(ter_val)
                self.stab_hard.append(hard_val)
                self.bin_life.append((bin_val0)/max(1,bin_val))
                self.ter_life.append((ter_val0)/max(1,ter_val))
                self.hard_life.append((hard_val0)/max(1,hard_val))
                self.dedt.append(np.mean(temp_dedt))
                self.dadt.append(np.mean(temp_dadt))

        self.bin_life = np.asarray(self.bin_life)
        self.ter_life = np.asarray(self.ter_life)
        self.hard_life = np.asarray(self.hard_life)
        self.stab_bin = np.asarray(self.stab_bin)
        self.stab_ter = np.asarray(self.stab_ter)
        self.stab_hard = np.asarray(self.stab_hard)
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
        hard_formed = np.empty(len(ini_pop))
        bin_life = np.empty(len(ini_pop))
        ter_life = np.empty(len(ini_pop))
        hard_life = np.empty(len(ini_pop))

        iter = -1
        for pop_ in ini_pop:
            iter += 1
            idx = np.where(self.pop == pop_)[0]
            bin_formed[iter] = np.mean(self.stab_bin[idx])
            ter_formed[iter] = np.mean(self.stab_ter[idx])
            hard_formed[iter] = np.mean(self.stab_hard[idx])
            bin_life[iter] = np.mean(self.bin_life[idx])
            ter_life[iter] = np.mean(self.ter_life[idx])
            hard_life[iter] = np.mean(self.hard_life[idx])
        bin_formed = np.asarray(bin_formed)
        ter_formed = np.asarray(ter_formed)
        hard_formed = np.asarray(hard_formed)
        bin_life = np.asarray(bin_life)
        ter_life = np.asarray(ter_life)
        hard_life = np.asarray(hard_life)

        norm_min = np.log10(min(np.min(bin_life[bin_formed>0]), np.min(ter_life[ter_formed>0]), np.min(hard_life[hard_formed>0])))
        norm_max = np.log10(max(np.max(bin_life[bin_formed>0]), np.max(ter_life[ter_formed>0]), np.max(hard_life[hard_formed>0])))

        fig, ax = plt.subplots()
        plot_ini.tickers_pop(ax, ini_pop)
        ax.set_xlabel(r'IMBH Population [$N$]')
        ax.set_ylabel(r'$N_{sys}$')
        normalise = plt.Normalize((norm_min), (norm_max))
        plt.scatter(ini_pop[bin_formed>0], bin_formed[bin_formed>0], edgecolors  = 'black', c = np.log10(bin_life[bin_formed>0]), norm = (normalise), label = 'Binary')
        plt.scatter(ini_pop[ter_formed>0], ter_formed[ter_formed>0], edgecolors  = 'black', c = np.log10(ter_life[ter_formed>0]), norm = (normalise), marker = 's', label = 'Triple')
        plt.scatter(ini_pop[hard_formed>0], hard_formed[hard_formed>0], edgecolors  = 'black', c = np.log10(hard_life[hard_formed>0]), norm = (normalise), marker = '^', label = 'Hard Binary')
        plt.colorbar().set_label(r'$\log_{10}\langle t_{sys} \rangle$ [kyr]')
        ax.set_xlim(5,105)
        ax.legend()
        plt.savefig('figures/sys_formation_N_plot.pdf', dpi=300, bbox_inches='tight')

        fig = plt.figure(figsize=(12.5, 6))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax1.set_xlabel(r'$\log_{10}\langle \dot{a} \rangle$ [pc Myr$^{-1}$]')
        ax1.set_ylabel(r'$\log_{10}\langle \dot{(1-e)} \rangle$ [Myr$^{-1}$]')
        ax2.set_xlabel(r'$\langle \dot{a} \rangle$ [pc Myr$^{-1}$]')
        ax2.set_ylabel(r'$\langle \dot{(1-e)} \rangle$ [Myr$^{-1}$]')
        for ax_ in [ax1, ax2]:
            plot_ini.tickers(ax)
            ax_.yaxis.set_major_formatter(mtick_formatter)
            ax_.xaxis.set_major_formatter(mtick_formatter)
        colour_axes = ax1.scatter(np.log10(10**3*self.dadt), np.log10(10**3*self.dedt), edgecolors = 'black', c = self.pop)
        colour_axes = ax2.scatter((10**3*self.dadt), (10**3*self.dedt), edgecolors = 'black', c = self.pop)
        plt.colorbar(colour_axes, ax = ax2, label = r'Initial Population')
        plt.savefig('figures/dadt_dedt_plot.pdf', dpi=300, bbox_inches='tight')

cst = sustainable_sys()
cst.system_formation()