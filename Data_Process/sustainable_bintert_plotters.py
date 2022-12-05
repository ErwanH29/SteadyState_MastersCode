from amuse.lab import *
from file_logistics import *
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import scipy.optimize
import warnings

class sustainable_sys(object):
    """
    Class which extracts ALL information of particle_trajectory files in an memory-optimised manner
    and appends the required data values into arrays.
    """

    def __init__(self):
        warnings.filterwarnings("ignore", category=RuntimeWarning) 

        filenameH = glob.glob(os.path.join('/media/erwanh/Elements/particle_trajectory_temp/*'))
        filenameGRX = glob.glob('/media/erwanh/Elements/GRX/particle_trajectory/*')
        filename = [natsort.natsorted(filenameH), natsort.natsorted(filenameGRX)]

        self.pop = [[ ], [ ]]
        self.dedt = [[ ], [ ]]
        self.dadt = [[ ], [ ]]
        self.semi_NN_avg = [[ ], [ ]]
        self.semi_NN_min = [[ ], [ ]]
        self.semi_t_avg = [[ ], [ ]]
        self.semi_t_min = [[ ], [ ]]
        self.sys_bin = [[ ], [ ]]
        self.sys_ter = [[ ], [ ]]
        self.bsys_time = [[ ], [ ]]
        self.tsys_time = [[ ], [ ]]
        self.bform_time = [[ ], [ ]]
        self.tform_time = [[ ], [ ]]
        self.pop_tracker = [[ ], [ ]]
        
        for int_ in range(2):
            for file_ in range(len(filename[int_])):
                with open(filename[int_][file_], 'rb') as input_file:
                    print('Reading file', file_, ':', input_file)
                    data = pkl.load(input_file)

                    pop = 10*round(0.1*np.shape(data)[0])
                    self.pop[int_].append(pop)

                    temp_dedt = [] 
                    temp_dadt = []
                    bsys_time = []
                    tsys_time = []
                    bin_key = [] 
                    ter_key = []

                    bin_val = 0 
                    ter_val = 0
                    bin_sys = 0 
                    ter_sys = 0
                    for parti_ in range(np.shape(data)[0]):
                        bin_sys = False
                        ter_sys = False
                        semi_nn_avg = []
                        semi_t_avg = []
                        self.bform_time[int_].append(-5)
                        self.tform_time[int_].append(-5)
                        self.pop_tracker[int_].append(10*round(0.1*np.shape(data)[0]))
                        if parti_ != 0:
                            temp_dedt.append((data.iloc[parti_][-2][8][0] - data.iloc[parti_][2][8][0])/(np.shape(data)[1]-3))
                            temp_dadt.append((data.iloc[parti_][-2][7][0]**-1 - data.iloc[parti_][2][7][0]**-1).value_in((units.pc)**-1)/(np.shape(data)[1]-3))
                            for col_ in range(np.shape(data)[1]-1):
                                nn_semi = data.iloc[parti_][col_][7][1]
                                if abs(data.iloc[parti_][col_][8][1]) < 1 \
                                    and nn_semi < 0.02 | units.parsec:
                                    mass1 = data.iloc[parti_][0][1]
                                    for part_ in range(np.shape(data)[0]):
                                        if data.iloc[part_][0][0] == data.iloc[parti_][col_][6][1]:
                                            mass2 = data.iloc[part_][0][1]
                                            bin_BE = ((constants.G*mass1*mass2)/(nn_semi)).value_in(units.J) 
                                            if bin_BE > (150000)**2*(1+mass1/mass2):   #Hard binary conditions based on Quinlan 1996b
                                                if not (bin_sys):                     #First formation time
                                                    formation_time = col_*1000
                                                    self.bform_time[int_][-1] = formation_time
                                                    bin_sys = True
                                                bin_val += 1
                                                bin_key.append(data.iloc[parti_][col_][6][1])
                                                bsys_time.append(col_)
                                                semi_nn_avg.append(nn_semi.value_in(units.pc))

                                                semi_outer = data.iloc[parti_][col_][7][2]
                                                ecc_outer = data.iloc[parti_][col_][8][2]
                                                #Calculate tertiary. The stability equality is based on Mardling and Aarseth 2001
                                                if ecc_outer < 1 and semi_outer < 0.1 | units.parsec:
                                                    for part_ in range(np.shape(data)[0]):
                                                        if data.iloc[part_][0][0] == data.iloc[parti_][col_][6][2]:
                                                            mass_outer = data.iloc[part_][0][1]

                                                    semi_ratio = semi_outer / nn_semi
                                                    equality = 2.8 * ((1+mass_outer/(mass1+mass2))*(1+ecc_outer)/(1-ecc_outer**2)**0.5)**0.4
                                                    if semi_ratio > equality:
                                                        if not (ter_sys):
                                                            formation_time = col_*1000
                                                            self.tform_time[int_][-1] = formation_time
                                                            ter_sys = True
                                                        ter_val += 1
                                                        ter_key.append([i for i in data.iloc[parti_][col_][6][2]][0])
                                                        tsys_time.append(col_)
                                                        semi_t_avg.append(semi_outer.value_in(units.pc))

                        if len(semi_nn_avg) > 0:
                            self.semi_NN_avg[int_].append(np.mean(semi_nn_avg))
                            self.semi_NN_min[int_].append(np.min(semi_nn_avg))
                        else:
                            self.semi_NN_avg[int_].append(-5)
                            self.semi_NN_min[int_].append(-5)
                        if len(semi_t_avg) > 0:
                            self.semi_t_avg[int_].append(np.mean(semi_t_avg))
                            self.semi_t_min[int_].append(np.min(semi_t_avg))
                        else:
                            self.semi_t_avg[int_].append(-5)
                            self.semi_t_min[int_].append(-5)

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
            self.semi_NN_avg[int_] = np.asarray([i for i in self.semi_NN_avg[int_]])
            self.semi_NN_min[int_] = np.asarray([i for i in self.semi_NN_min[int_]])
            self.semi_t_avg[int_] = np.asarray([i for i in self.semi_t_avg[int_]])
            self.semi_t_min[int_] = np.asarray([i for i in self.semi_t_min[int_]])
            self.pop[int_] = np.asarray([i for i in self.pop[int_]])
            self.bform_time[int_] = np.asarray([i for i in self.bform_time[int_]])
            self.tform_time[int_] = np.asarray([i for i in self.tform_time[int_]])

        self.semi_NN_max = []
        self.semi_t_max = []
        with open('figures/binary_hierarchical/output/system_summary.txt', 'w') as file:
            integrator = ['Hermite', 'GRX']
            for int_ in range(2):
                semi_NN_avg = [ ]
                semi_NN_min = [ ]
                semi_t_avg = [ ]
                semi_t_min = [ ]
                pop_arr = np.unique(self.pop_tracker[int_])
                for pop_ in pop_arr:
                    idx = np.argwhere(self.pop_tracker[int_] == pop_)
                    semi_av = self.semi_NN_avg[int_][idx][self.semi_NN_avg[int_][idx] > 0]
                    semi_mv = self.semi_NN_min[int_][idx][self.semi_NN_min[int_][idx] > 0]
                    semi_avt = self.semi_t_avg[int_][idx][self.semi_t_avg[int_][idx] > 0]
                    semi_mvt = self.semi_t_min[int_][idx][self.semi_t_min[int_][idx] > 0]
                    self.semi_NN_max.append(np.mean(semi_av))
                    self.semi_t_max.append(np.mean(semi_avt))
                    if len(semi_av) > 0:
                        semi_NN_avg.append(np.mean(semi_av))
                        semi_NN_min.append(np.min(semi_mv))
                    if len(semi_avt) > 0:
                        semi_t_avg.append(np.mean(semi_avt))
                        semi_t_min.append(np.min(semi_mvt))
                file.write('\nData for '+str(integrator[int_]))
                file.write('\nAverage binary semi-major axis per population:   '+str(pop_arr)+' : '+str(semi_NN_avg))
                file.write('\nMinimum binary semi-major axis per population:   '+str(pop_arr)+' : '+str(semi_NN_min))
                file.write('\nAverage tertiary semi-major axis per population: '+str(pop_arr)+' : '+str(semi_t_avg))
                file.write('\nMinimum tertiary semi-major axis per population: '+str(pop_arr)+' : '+str(semi_t_min))

    def system_formation_data(self, int_, ini_pop):
        """
        Extract data to form plots
        
        Inputs:
        int_:    Integrator defining which data is being used
        ini_pop: The population of the given simulation results
        """

        bin_formed = np.empty(len(ini_pop))
        ter_formed = np.empty(len(ini_pop))
        bsys_time = np.empty(len(ini_pop))
        tsys_time = np.empty(len(ini_pop))
        bform_time = np.empty(len(ini_pop))
        tform_time = np.empty(len(ini_pop))
        semi_major_bin = np.empty(len(ini_pop))
        semi_major_ter = np.empty(len(ini_pop))

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
        
            idx3 = np.argwhere(self.pop_tracker[int_] == pop_)
            semi_av = self.semi_NN_avg[int_][idx3][self.semi_NN_avg[int_][idx3] > 0]
            semi_avt = self.semi_t_avg[int_][idx3][self.semi_t_avg[int_][idx3] > 0]
            if len(semi_av) > 0:
                semi_major_bin[iter] = np.mean(semi_av)
            if len(semi_avt) > 0:
                semi_major_ter[iter] = np.mean(semi_avt)

        with open('figures/binary_hierarchical/output/'+str(integrator[int_])+'bin_ter_systems.txt', 'w') as file:
            file.write(str(integrator[int_])+' first binary avg. formation time')
            for pop_ in range(len(bform_time)):
                file.write('\nPopulation: '+str(ini_pop[pop_])+': '+str(bform_time[pop_]/10**3)+' kyr')
            file.write('\n\nfirst tertiary avg. formation time')
            for pop_ in range(len(tform_time)):
                file.write('\nPopulation: '+str(ini_pop[pop_])+': '+str(tform_time[pop_]/10**3)+' kyr')

        return bin_formed, ter_formed, bsys_time, tsys_time, semi_major_bin, semi_major_bin
    
    def system_formation_plotter(self):
        """
        Function to plot various 'sustainable system' plots
        """

        def log_fit(xval, slope, alpha, yint):
            return (slope)**(xval+alpha) + yint

        plot_ini = plotter_setup()
        mtick_formatter = mtick.FormatStrFormatter('%0.2f')
        integrator = ['Hermite', 'GRX']

        p21ymin = min(np.nanmin(np.log10(10**6*(self.dedt[0]))), np.nanmin(np.log10(10**6*(self.dedt[1]))))
        p21ymax = max(np.nanmax(np.log10(10**6*(self.dedt[0]))), np.nanmax(np.log10(10**6*(self.dedt[1]))))
        p21xmin = min(np.nanmin(np.log10(10**6*(self.dadt[0]))), np.nanmin(np.log10(10**6*(self.dadt[1]))))
        p21xmax = max(np.nanmax(np.log10(10**6*(self.dadt[0]))), np.nanmax(np.log10(10**6*(self.dadt[1]))))
        normalise_p1 = plt.Normalize(0, max(max(self.semi_NN_max), max(self.semi_t_max)))
        normalise_p2 = plt.Normalize(10, 100)
    
        fig = plt.figure(figsize=(14, 6))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax_ = [ax1, ax2]
        for int_ in range(1):
            ini_pop = np.unique(self.pop[int_])
            plot_ini.tickers_pop(ax_[int_], ini_pop)
            ax_[int_].set_title(integrator[int_])
            ax_[int_].set_xlabel(r'IMBH Population [$N$]')
            ax_[int_].set_ylabel(r'$\log_{10} t_{\rm{sys}} / t_{\rm{sim}}$')
            ax_[int_].set_ylim(-5.2, 0)
            bin_formed, ter_formed, bsys_time, tsys_time, semi_bin, semi_ter = self.system_formation_data(int_, ini_pop)
            colour_axes = ax_[int_].scatter(ini_pop[bin_formed>0], np.log10(bsys_time[bin_formed>0]), edgecolors  = 'black', c = (semi_bin[bin_formed>0]), norm = (normalise_p1), label = 'Stable Binary')
            ax_[int_].scatter(ini_pop[ter_formed>0], np.log10(tsys_time[ter_formed>0]), edgecolors  = 'black', c = (semi_ter[ter_formed>0]), norm = (normalise_p1), marker = 's', label = 'Stable Triple')

            """p0 = (10, 2, 0)
            params_bin, cv = scipy.optimize.curve_fit(log_fit, ini_pop[bin_formed>0], bsys_time[bin_formed>0], p0, maxfev = 2000)
            slope_bin, alpha_bin, intercept_bin = params_bin
            params_ter, cv = scipy.optimize.curve_fit(log_fit, ini_pop[ter_formed>0], bsys_time[ter_formed>0], p0, maxfev = 2000)
        
        y_arr = ([log_fit(i, slope_bin, alpha_bin, intercept_bin) for i in x_arr])
        print(slope_bin, alpha_bin, intercept_bin)"""

        p0 = (1, 1, 10**-5)
        params, cv = scipy.optimize.curve_fit(log_fit, ini_pop[bin_formed>0], (bsys_time[bin_formed>0]), p0, maxfev = 10000, method = 'trf')
        slope, beta, yint = params

        x_arr = np.linspace(10,100)
        pb = np.polyfit(ini_pop[bin_formed>0], (bsys_time[bin_formed>0]), 2)
        fit_pb = np.polynomial.polynomial.Polynomial(pb[::-1])(x_arr)
        pt = np.polyfit(ini_pop[ter_formed>0], (tsys_time[ter_formed>0]), 2)
        fit_pt = np.polynomial.polynomial.Polynomial(pt[::-1])(x_arr)
        ax1.plot(x_arr, np.log10(fit_pb))
        ax1.plot(x_arr, np.log10(fit_pt))
        plt.colorbar(colour_axes, ax=ax2, label = r'$\langle N_{\rm{sys}} \rangle$ ')
        ax2.legend()
        plt.savefig('figures/binary_hierarchical/sys_formation_N_plot.pdf', dpi=300, bbox_inches='tight')

        print('Slope : ', pb[0])
        print('y-int : ', pb[1])

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
            colour_axes = ax1.scatter(np.log10(10**6*self.dadt[int_]), np.log10(10**6*self.dedt[int_]),
                                      norm = normalise_p2, edgecolors = 'black', c = self.pop[int_])
            colour_axes = ax2.scatter(10**6*self.dadt[int_], (10**6*self.dedt[int_]), 
                                      norm = normalise_p2, edgecolors = 'black', c = self.pop[int_])

            plt.colorbar(colour_axes, ax = ax2, label = r'Initial Population')
            plt.savefig('figures/binary_hierarchical/simavg_dadt_dedt_plot'+str(integrator[int_])+'.pdf', dpi=300, bbox_inches='tight')
            plt.clf()
