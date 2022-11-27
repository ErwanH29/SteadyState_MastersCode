from amuse.lab import *
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
        self.sys_bin = [[ ], [ ]]
        self.sys_ter = [[ ], [ ]]
        self.bform_time = [[ ], [ ]]
        self.tform_time = [[ ], [ ]]
        self.bsys_time = [[ ], [ ]]
        self.tsys_time = [[ ], [ ]]
        self.pop_tracker = [[ ], [ ]]
        
        for int_ in range(2):
            for file_ in range(len(filename[int_])):
                with open(filename[int_][file_], 'rb') as input_file:
                    print('Reading file', file_, ':', input_file)
                    data = pkl.load(input_file)

                    pop = np.shape(data)[0]
                    if pop %10 != 0:
                        pop -=1
                    self.pop[int_].append(pop)
                    bin_key = [] ; ter_key = []
                    temp_dedt = [] ; temp_dadt = []
                    bin_val = 0 ; ter_val = 0
                    bin_sys = 0 ; ter_sys = 0
                    bsys_time = []
                    tsys_time = []

                    avg_mass = 0
                    for parti_ in range(np.shape(data)[0]-1):
                        avg_mass +=  data.iloc[parti_+1][0][1].value_in(units.kg)
                    avg_mass /= (np.shape(data)[0])

                    for parti_ in range(np.shape(data)[0]):
                        bin_sys = False
                        ter_sys = False
                        self.bform_time[int_].append(-5)
                        self.tform_time[int_].append(-5)
                        self.pop_tracker[int_].append(10*round(0.1*np.shape(data)[0]))
                        if parti_ != 0:
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
                                                if not (bin_sys):                     #First formation time
                                                    formation_time = col_*1000
                                                    self.bform_time[int_][-1] = formation_time
                                                    bin_sys = True
                                                bin_val += 1
                                                bin_key.append([i for i in data.iloc[parti_][col_][6][1]][0])
                                                bsys_time.append(col_)

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
                                                        if not (ter_sys):
                                                            formation_time = col_*1000
                                                            self.tform_time[int_][-1] = formation_time
                                                            ter_sys = True

                                                        ter_val += 1
                                                        ter_key.append([i for i in data.iloc[parti_][col_][6][2]][0])
                                                        tsys_time.append(col_)

                        if parti_ == (np.shape(data)[0]-1):
                            if len(bin_key) > 0:
                                bin_formed = len(np.unique(bin_key))
                                bin_sys += bin_formed
                            if len(ter_key) > 0:
                                ter_formed = len(np.unique(ter_key))
                                ter_sys += ter_formed
                                
                    self.bsys_time[int_].append((len(np.unique(bsys_time)))/(col_))
                    self.tsys_time[int_].append((len(np.unique(tsys_time)))/(col_))
                    self.sys_bin[int_].append(bin_sys)
                    self.sys_ter[int_].append(ter_sys)
                    self.dedt[int_].append(np.mean(temp_dedt))
                    self.dadt[int_].append(np.mean(temp_dadt))

        for int_ in range(2):
            self.bsys_time[int_] = np.asarray([i for i in self.bsys_time[int_]])
            self.tsys_time[int_] = np.asarray([i for i in self.tsys_time[int_]])
            self.sys_bin[int_] = np.asarray([i for i in self.sys_bin[int_]])
            self.sys_ter[int_] = np.asarray([i for i in self.sys_ter[int_]])
            self.dedt[int_] = np.asarray([i for i in self.dedt[int_]])
            self.dadt[int_] = np.asarray([i for i in self.dadt[int_]])
            self.pop[int_] = np.asarray([i for i in self.pop[int_]])
            self.bform_time[int_] = np.asarray([i for i in self.bform_time[int_]])
            self.tform_time[int_] = np.asarray([i for i in self.tform_time[int_]])

    def system_formation_data(self, int_, ini_pop):
        """
        Function to plot various 'sustainable system' plots
        """

        bin_formed = np.empty(len(ini_pop))
        ter_formed = np.empty(len(ini_pop))
        bsys_time = np.empty(len(ini_pop))
        tsys_time = np.empty(len(ini_pop))
        bform_time = np.empty(len(ini_pop))
        tform_time = np.empty(len(ini_pop))

        integrator = ['Hermite', 'GRX']
        iter = -1
        for pop_ in ini_pop:
            iter += 1
            idx = np.where(self.pop[int_] == pop_)[0]
            idx2 = np.where(self.pop_tracker[int_] == pop_)[0]
            
            bin_formed[iter] = np.mean(self.sys_bin[int_][idx])
            ter_formed[iter] = np.mean(self.sys_ter[int_][idx])
            bsys_time[iter] = np.mean(self.bsys_time[int_][idx])
            tsys_time[iter] = np.mean(self.tsys_time[int_][idx])
            bform_time[iter] = np.mean(np.asarray(self.bform_time[int_][idx2])[(self.bform_time[int_][idx2]) >= 0])
            tform_time[iter] = np.mean(np.asarray(self.tform_time[int_][idx2])[(self.tform_time[int_][idx2]) >= 0])

        with open('figures/binary_hierarchical/output/'+str(integrator[int_])+'bin_ter_systems.txt', 'w') as file:
            file.write(str(integrator[int_])+' first binary avg. formation time')
            for pop_ in range(len(bform_time)):
                file.write('\nPopulation: '+str(ini_pop[pop_])+': '+str(bform_time[pop_]/10**3)+' kyr')
            file.write('\n\nfirst tertiary avg. formation time')
            for pop_ in range(len(tform_time)):
                file.write('\nPopulation: '+str(ini_pop[pop_])+': '+str(tform_time[pop_]/10**3)+' kyr')

        return bin_formed, ter_formed, bsys_time, tsys_time
    
    def system_formation_plotter(self):

        plot_ini = plotter_setup()
        mtick_formatter = mtick.FormatStrFormatter('%0.2f')
        integrator = ['Hermite', 'GRX']
        norm_min = min(np.nanmin(np.mean(self.sys_bin[0])), np.nanmin(np.mean(self.sys_bin[1])),
                       np.nanmin(np.mean(self.sys_ter[0])), np.nanmin(np.mean(self.sys_ter[1])))
        norm_max = max(np.nanmax(np.mean(self.sys_bin[0])), np.nanmax(np.mean(self.sys_bin[1])),
                       np.nanmax(np.mean(self.sys_ter[0])), np.nanmin(np.mean(self.sys_ter[1])))
        normalise_p1 = plt.Normalize((norm_min), (norm_max))
        normalise_p2 = plt.Normalize(10, 100)

        p21ymin = min(np.nanmin(np.log10(10**6*(self.dedt[0]))), np.nanmin(np.log10(10**6*(self.dedt[1]))))
        p21ymax = max(np.nanmax(np.log10(10**6*(self.dedt[0]))), np.nanmax(np.log10(10**6*(self.dedt[1]))))
        p21xmin = min(np.nanmin(np.log10(10**6*(self.dadt[0]))), np.nanmin(np.log10(10**6*(self.dadt[1]))))
        p21xmax = max(np.nanmax(np.log10(10**6*(self.dadt[0]))), np.nanmax(np.log10(10**6*(self.dadt[1]))))
    
        fig = plt.figure(figsize=(14, 6))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax_ = [ax1, ax2]
        for int_ in range(2):
            ini_pop = np.unique(self.pop[int_])
            plot_ini.tickers_pop(ax_[int_], ini_pop)
            ax_[int_].set_title(integrator[int_])
            ax_[int_].set_xlabel(r'IMBH Population [$N$]')
            ax_[int_].set_ylabel(r'$t_{\rm{sys}} / t_{\rm{sim}}$')
            ax_[int_].set_ylim(0, 1.03)
            bin_formed, ter_formed, bsys_time, tsys_time = self.system_formation_data(int_, ini_pop)
            colour_axes = ax_[int_].scatter(ini_pop[bin_formed>0], np.log10(bsys_time[bin_formed>0]), edgecolors  = 'black', c = (bin_formed[bin_formed>0]), norm = (normalise_p1), label = 'Stable Binary')
            ax_[int_].scatter(ini_pop[ter_formed>0], np.log10(tsys_time[ter_formed>0]), edgecolors  = 'black', c = (ter_formed[ter_formed>0]), norm = (normalise_p1), marker = 's', label = 'Stable Triple')
        plt.colorbar(colour_axes, ax=ax2, rotation=270, label = r'$\log_{10}\langle N_{\rm{sys}} \rangle$ ')
        ax2.legend()
        plt.savefig('figures/binary_hierarchical/sys_formation_N_plot.pdf', dpi=300, bbox_inches='tight')

        for int_ in range(2):
            fig = plt.figure(figsize=(12.5, 6))
            ax1 = fig.add_subplot(121)
            ax2 = fig.add_subplot(122)
            for ax_ in [ax1, ax2]:
                plot_ini.tickers(ax_, 'plot')
            ax1.set_title(integrator[int_])
            ax1.yaxis.set_major_formatter(mtick_formatter)            
            ax1.set_xlim(0.9*p21xmin, 1.1*p21xmax)
            ax1.set_ylim(0.9*p21ymin, 1.1*p21ymax)        
            ax1.set_ylabel(r'$\log_{10}\langle \dot{(1-e)} \rangle_{SMBH}$ [Gyr$^{-1}$]')
            ax2.set_ylabel(r'$\langle \dot{(1-e)} \rangle_{SMBH}$ [Gyr$^{-1}$]')
            ax1.set_xlabel(r'$\log_{10}\langle \dot{a}^{-1} \rangle_{SMBH}$ [pc$^{-1}$Gyr$^{-1}$]')
            ax2.set_xlabel(r'$\langle \dot{a}^{-1} \rangle_{SMBH}$ [pc$^{-1}$Gyr$^{-1}$]')
            ax1.xaxis.set_major_formatter(mtick_formatter)
            colour_axes = ax1.scatter(np.log10(10**6*self.dadt[int_]), np.log10(10**6*self.dedt[int_]), cmap='tab10', 
                                      norm = normalise_p2, edgecolors = 'black', c = self.pop[int_])
            colour_axes = ax2.scatter(10**6*self.dadt[int_], (10**6*self.dedt[int_]), cmap='tab10', 
                                      norm = normalise_p2, edgecolors = 'black', c = self.pop[int_])

            plt.colorbar(colour_axes, ax = ax2, rotation=270, label = r'Initial Population')
            plt.savefig('figures/binary_hierarchical/dadt_dedt_plot'+str(integrator[int_])+'.pdf', dpi=300, bbox_inches='tight')
            plt.clf()
