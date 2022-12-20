from amuse.lab import *
from file_logistics import *
from tGW_plotters import *
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

        GW_calcs = gw_calcs()

        filenameH = glob.glob(os.path.join('/media/erwanh/Elements/Hermite/particle_trajectory_temp/*'))
        filenameGRX = glob.glob('/media/erwanh/Elements/GRX/particle_trajectory_temp/*')
        filename = [natsort.natsorted(filenameH), natsort.natsorted(filenameGRX)]

        self.pop = [[ ], [ ]]
        self.pop_bin = [[ ], [ ]]
        self.pop_ter = [[ ], [ ]]
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
        self.parti_velb = [[ ], [ ]]
        self.GW_freqbin = [[ ], [ ]]
        self.GW_strainbin = [[ ], [ ]]
        self.GW_timeb = [[ ], [ ]]
        self.parti_velt = [[ ], [ ]]
        self.GW_freqter = [[ ], [ ]]
        self.GW_strainter = [[ ], [ ]]
        self.GW_timet = [[ ], [ ]]
        self.GW_bmass = [[ ], [ ]]
        self.GW_tmass = [[ ], [ ]]
        self.hard_bin = [[ ], [ ]]
        self.hard_ter = [[ ], [ ]]
        
        for int_ in range(2):
            for file_ in range(len(filename[int_])):
                with open(filename[int_][file_], 'rb') as input_file:
                    print('Reading file', file_, ':', input_file)
                    file_size = os.path.getsize(filename[int_][file_])
                    if file_size < 2.8e9:
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
                        if pop > 5 and pop <= 50:
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
                                        nn_ecc = data.iloc[parti_][col_][8][1]
                                        mass1 = data.iloc[parti_][0][1]
                                        for part_ in range(np.shape(data)[0]):
                                            if data.iloc[part_][0][0] == data.iloc[parti_][col_][6][1]:
                                                mass2 = data.iloc[part_][0][1]
                                                mass = max(mass1, mass2)
                                                hard = False
                                                bin = False
                                                print((constants.G*mass)/(4*(150000*(1 | units.ms))**2))
                                                if  nn_semi < (constants.G*mass)/(4*(150000*(1 | units.ms))**2):   #Hard binary conditions based on Quinlan 1996b
                                                    hard = True
                                                    bin = True
                                                    self.hard_bin[int_].append(1)
                                                print(((constants.G*mass)/(4*(15000*(1 | units.ms))**2)).value_in(units.pc))
                                                STOP
                                                if not (hard) and nn_semi < (constants.G*mass)/(4*(1500*(1 | units.ms))**2) and abs(nn_ecc) < 1:  #Value chosen for 1000AU
                                                    self.hard_bin[int_].append(-5)
                                                    bin = True

                                                if (bin):
                                                    if not (bin_sys):  #First formation time
                                                        formation_time = col_*1000
                                                        self.bform_time[int_][-1] = formation_time
                                                        bin_sys = True
                                                    bin_val += 1
                                                    bin_key.append(data.iloc[parti_][col_][6][1])
                                                    bsys_time.append(col_)
                                                    semi_nn_avg.append(nn_semi.value_in(units.pc))
                                                    self.pop_bin[int_].append(pop)

                                                    strain = GW_calcs.gw_strain(nn_semi, nn_ecc, mass1, mass2)
                                                    freq = GW_calcs.gw_freq(nn_semi, nn_ecc, mass1, mass2)
                                                    GW_time = GW_calcs.gw_timescale(nn_semi, nn_ecc, mass1, mass2)
                                                    self.GW_strainbin[int_].append(strain)
                                                    self.GW_freqbin[int_].append(freq)
                                                    self.GW_timeb[int_].append(float(GW_time.value_in(units.Myr)))

                                                    velx = data.iloc[parti_][col_][3][0]
                                                    vely = data.iloc[parti_][col_][3][1]
                                                    velz = data.iloc[parti_][col_][3][2]
                                                    vel = np.sqrt(velx**2+vely**2+velz**2)
                                                    self.parti_velb[int_].append(float(vel.value_in(units.kms)))
                                                    self.GW_bmass[int_].append(mass2.value_in(units.MSun))

                                                    semi_outer = data.iloc[parti_][col_][7][2]
                                                    ecc_outer = data.iloc[parti_][col_][8][2]

                                                    #Calculate tertiary. The stability equality is based on Mardling and Aarseth 2001
                                                    if ecc_outer < 1 and semi_outer < 0.1 | units.parsec:
                                                        for part_ in range(np.shape(data)[0]):
                                                            if data.iloc[part_][0][0] == data.iloc[parti_][col_][6][2]:
                                                                mass_outer = data.iloc[part_][0][1]

                                                        semi_ratio = semi_outer / nn_semi
                                                        equality = 2.8 * ((1+mass_outer/(mass1+mass2))*(1+ecc_outer)/(1-ecc_outer)**0.5)**0.4
                                                        if semi_ratio > equality:
                                                            if not (ter_sys):
                                                                formation_time = col_*1000
                                                                self.tform_time[int_][-1] = formation_time
                                                                ter_sys = True
                                                            ter_val += 1
                                                            ter_key.append([i for i in data.iloc[parti_][col_][6][2]][0])
                                                            tsys_time.append(col_)
                                                            semi_t_avg.append(semi_outer.value_in(units.pc))
                                                            self.pop_ter[int_].append(pop)

                                                            GW_time = GW_calcs.gw_timescale(semi_outer, ecc_outer, mass1, mass_outer)

                                                            strain = GW_calcs.gw_strain(semi_outer, ecc_outer, mass1, mass_outer)
                                                            freq = GW_calcs.gw_freq(semi_outer, ecc_outer, mass1, mass_outer)
                                                            self.GW_strainter[int_].append(strain)
                                                            self.GW_freqter[int_].append(freq)
                                                            self.parti_velt[int_].append(float(vel.value_in(units.kms)))
                                                            self.GW_timet[int_].append(float(GW_time.value_in(units.Myr)))
                                                            self.GW_tmass[int_].append([mass2.value_in(units.MSun), mass_outer.value_in(units.MSun)])

                                                            if (hard):
                                                                self.hard_ter[int_].append(1)
                                                            else:
                                                                self.hard_ter[int_].append(-5)

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
            self.pop_bin[int_] = np.asarray(self.pop_bin[int_])
            self.pop_ter[int_] = np.asarray(self.pop_ter[int_])
            self.bsys_time[int_] = np.asarray(self.bsys_time[int_])
            self.tsys_time[int_] = np.asarray(self.tsys_time[int_])
            self.sys_bin[int_] = np.asarray(self.sys_bin[int_])
            self.sys_ter[int_] = np.asarray(self.sys_ter[int_])
            self.dedt[int_] = np.asarray(self.dedt[int_])
            self.dadt[int_] = np.asarray(self.dadt[int_])
            self.semi_NN_avg[int_] = np.asarray(self.semi_NN_avg[int_])
            self.semi_NN_min[int_] = np.asarray(self.semi_NN_min[int_])
            self.semi_t_avg[int_] = np.asarray(self.semi_t_avg[int_])
            self.semi_t_min[int_] = np.asarray(self.semi_t_min[int_])
            self.pop[int_] = np.asarray(self.pop[int_])
            self.bform_time[int_] = np.asarray(self.bform_time[int_])
            self.tform_time[int_] = np.asarray(self.tform_time[int_])
            self.parti_velb[int_] = np.asarray(self.parti_velb[int_])
            self.GW_freqbin[int_] = np.asarray(self.GW_freqbin[int_])
            self.GW_strainbin[int_] = np.asarray(self.GW_strainbin[int_])
            self.GW_timeb[int_] = np.asarray(self.GW_timeb[int_])
            self.GW_bmass[int_] = np.asarray(self.GW_bmass[int_])
            self.parti_velt[int_] = np.asarray(self.parti_velt[int_])
            self.GW_freqter[int_] = np.asarray(self.GW_freqter[int_])
            self.GW_strainter[int_] = np.asarray(self.GW_strainter[int_])
            self.GW_timet[int_] = np.asarray(self.GW_timet[int_])
            self.GW_tmass[int_] = np.asarray(self.GW_tmass[int_])
            self.hard_ter[int_] = np.asarray(self.hard_ter[int_])
            self.hard_bin[int_] = np.asarray(self.hard_bin[int_])

        self.semi_NN_max = []
        self.semi_t_max = []
        with open('figures/binary_hierarchical/output/system_summary.txt', 'w') as file:
            self.integrator = ['Hermite', 'GRX']
            for int_ in range(2):
                semi_NN_avg = [ ]
                semi_NN_min = [ ]
                semi_t_avg = [ ]
                semi_t_min = [ ]
                pop_arr = np.unique(self.pop_tracker[int_])
                file.write('\n\nData for '+str(self.integrator[int_]+' in pc'))
                for pop_ in pop_arr:
                    idx = np.argwhere(self.pop_tracker[int_] == pop_)
                    semi_av = self.semi_NN_avg[int_][idx][self.semi_NN_avg[int_][idx] > 0]
                    semi_mv = self.semi_NN_min[int_][idx][self.semi_NN_min[int_][idx] > 0]
                    semi_avt = self.semi_t_avg[int_][idx][self.semi_t_avg[int_][idx] > 0]
                    semi_mvt = self.semi_t_min[int_][idx][self.semi_t_min[int_][idx] > 0]
                    self.semi_NN_max.append('{:.7f}'.format(np.mean(semi_av)))
                    self.semi_t_max.append('{:.7f}'.format(np.mean(semi_avt)))
                    if len(semi_av) > 0:
                        semi_NN_avg.append('{:.7f}'.format(np.mean(semi_av)))
                        semi_NN_min.append('{:.7f}'.format(np.min(semi_mv)))
                    if len(semi_avt) > 0:
                        semi_t_avg.append('{:.7f}'.format(np.mean(semi_avt)))
                        semi_t_min.append('{:.7f}'.format(np.min(semi_mvt)))
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

        GW_calcs = gw_calcs()
        tH = GW_calcs.tH

        bin_formed = np.empty(len(ini_pop))
        ter_formed = np.empty(len(ini_pop))
        bsys_time = np.empty(len(ini_pop))
        tsys_time = np.empty(len(ini_pop))
        bform_time = np.empty(len(ini_pop))
        bform_ini = np.empty(len(ini_pop))
        tform_time = np.empty(len(ini_pop))
        tform_ini = np.empty(len(ini_pop))
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
            bform_ini[iter] = len(np.asarray(self.bform_time[int_][idx2])[(self.bform_time[int_][idx2]) == 0])
            tform_time[iter] = np.mean(np.asarray(self.tform_time[int_][idx2])[(self.tform_time[int_][idx2]) >= 0])
            tform_ini[iter] = len(np.asarray(self.tform_time[int_][idx2])[(self.tform_time[int_][idx2]) == 0])
        
            semi_av = self.semi_NN_avg[int_][idx2][self.semi_NN_avg[int_][idx2] > 0]
            semi_avt = self.semi_t_avg[int_][idx2][self.semi_t_avg[int_][idx2] > 0]
            if len(semi_av) > 0:
                semi_major_bin[iter] = np.mean(semi_av)
            if len(semi_avt) > 0:
                semi_major_ter[iter] = np.mean(semi_avt)

        with open('figures/binary_hierarchical/output/'+str(integrator[int_])+'bin_ter_systems.txt', 'w') as file:
            file.write(str(integrator[int_])+' first binary avg. formation time')
            for pop_ in range(len(bform_time)):
                file.write('\nPopulation: '+str(ini_pop[pop_])+':         '+str(bform_time[pop_]/10**3)+' kyr')
            file.write('\n# Binary systems at initialisation:             '+str(bform_ini[pop_]))
            file.write('\n5 shortest binary GW merger time_scales [Myr]:  ')
            for i in range(5):
                file.write('{:.2f}'.format(np.sort(self.GW_timeb[int_][self.GW_timeb[int_] > 0])[i])+',   ')
            file.write('\nFraction of mergers within Hubble time:         '+str(len(self.GW_timeb[int_][self.GW_timeb[int_] < tH.value_in(units.Myr)])) +'/'+str(len(self.GW_timeb[int_])))
            file.write('\nFraction of IMBH-IMBH binaries:                 '+str(len(self.GW_bmass[int_][self.GW_bmass[int_] < 10**6])) +'/'+str(len(self.GW_bmass[int_])))

            file.write('\n\nfirst tertiary avg. formation time')
            for pop_ in range(len(tform_time)):
                file.write('\nPopulation: '+str(ini_pop[pop_])+':          '+str(tform_time[pop_]/10**3)+' kyr')
            file.write('\n# Tertiary systems at initialisation:            '+str(tform_ini[pop_]))
            file.write('\n5 shortest tertiary GW merger time scales [Myr]: ')
            for i in range(5):
                file.write('{:.2f}'.format(np.sort(self.GW_timet[int_][self.GW_timet[int_] > 0])[i])+',   ')
            file.write('\nFraction of mergers within Hubble time:          '+str(len(self.GW_timet[int_][self.GW_timet[int_] < tH.value_in(units.Myr)])) +'/'+str(len(self.GW_timet[int_])))
            file.write('\nFraction of (IMBH-IMBH)-SMBH tertiaries:         '+str(len(self.GW_tmass[int_][0][(self.GW_tmass[int_][0] < 10**6) & (self.GW_tmass[int_][1] > 10**6)])) +'/'+str(np.shape(self.GW_tmass[int_])[1]))
            file.write('\nFraction of (IMBH-SMBH)-IMBH tertiaries:         '+str(len(self.GW_tmass[int_][0][(self.GW_tmass[int_][1] < 10**6) & (self.GW_tmass[int_][0] > 10**6)])) +'/'+str(np.shape(self.GW_tmass[int_])[1]))
            file.write('\nFraction of (IMBH-IMBH)-IMBH tertiaries:         '+str(len(self.GW_tmass[int_][0][(self.GW_tmass[int_][1] < 10**6) & (self.GW_tmass[int_][0] < 10**6)])) +'/'+str(np.shape(self.GW_tmass[int_])[1]))

        return bin_formed, ter_formed, bsys_time, tsys_time, semi_major_bin, semi_major_bin
    
    def system_formation_plotter(self):
        """
        Function to plot various 'sustainable system' plots
        """

        def log_fit(xval, slope, alpha, yint):
            """
            Function of the form slope*10**(x+alpha)+yint
            """
            return slope*10**(xval+alpha) + yint

        plot_ini = plotter_setup()
        mtick_formatter = mtick.FormatStrFormatter('%0.2f')
        integrator = ['Hermite', 'GRX']

        """p21ymin = min(np.nanmin(np.log10(10**6*(self.dedt[0]))), np.nanmin(np.log10(10**6*(self.dedt[1]))))
        p21ymax = max(np.nanmax(np.log10(10**6*(self.dedt[0]))), np.nanmax(np.log10(10**6*(self.dedt[1]))))
        p21xmin = min(np.nanmin(np.log10(10**6*(self.dadt[0]))), np.nanmin(np.log10(10**6*(self.dadt[1]))))
        p21xmax = max(np.nanmax(np.log10(10**6*(self.dadt[0]))), np.nanmax(np.log10(10**6*(self.dadt[1]))))"""
        p21ymin = (np.nanmin(np.log10(10**6*(self.dedt[0]))))
        p21ymax = (np.nanmax(np.log10(10**6*(self.dedt[0]))))
        p21xmin = (np.nanmin(np.log10(10**6*(self.dadt[0]))))
        p21xmax = (np.nanmax(np.log10(10**6*(self.dadt[0]))))
        normalise_p1 = plt.Normalize(0, 500)
        normalise_p2 = plt.Normalize(10, 100)
    
        fig = plt.figure(figsize=(14, 6))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax_ = [ax1, ax2]
        for int_ in range(1):
            ini_pop = np.unique(self.pop[int_])
            ax_[int_].set_title(integrator[int_])
            ax_[int_].set_xlabel(r'IMBH Population [$N$]')
            ax_[int_].set_ylabel(r'$\log_{10} t_{\rm{sys}} / t_{\rm{sim}}$')
            ax_[int_].set_ylim(-5.2, 0)
            bin_formed, ter_formed, bsys_time, tsys_time, semi_bin, semi_ter = self.system_formation_data(int_, ini_pop)
            colour_axes = ax_[int_].scatter(ini_pop[bin_formed>0], np.log10(bsys_time[bin_formed>0]), edgecolors  = 'black', c = (bin_formed), norm = (normalise_p1), label = 'Stable Binary')
            ax_[int_].scatter(ini_pop[ter_formed>0], np.log10(tsys_time[ter_formed>0]), edgecolors  = 'black', c = (ter_formed), norm = (normalise_p1), marker = 's', label = 'Stable Triple')

            """p0 = (10, 2, 0)
            params_bin, cv = scipy.optimize.curve_fit(log_fit, ini_pop[bin_formed>0], bsys_time[bin_formed>0], p0, maxfev = 2000)
            slope_bin, alpha_bin, intercept_bin = params_bin
            params_ter, cv = scipy.optimize.curve_fit(log_fit, ini_pop[ter_formed>0], bsys_time[ter_formed>0], p0, maxfev = 2000)
            sloper_ter, alpha_ter, intercept_ter = params_ter
            #ax1.plot(x_arr, [np.log10(log_fit(i, slope_bin, alpha_bin, intercept_bin)) for i in x_arr])
            #ax1.plot(x_arr, [np.log10(log_fit(i, slope_ter, alpha_ter, intercept_ter)) for i in x_arr])
        print(slope_bin, alpha_bin, intercept_bin)"""

        plot_ini.tickers_pop(ax1, self.pop[0], 'Hermite')
        plot_ini.tickers_pop(ax2, self.pop[1], 'GRX')

        x_arr = np.linspace(10,100)
        #ax1.text(-5, 80, r'$N = \frac{1}{{{}}}(\log_{10}(t_{\rm{sys}} / t_{\rm{sim}})-{{}}-{{}}'.format(slope, yint, beta))
        plt.colorbar(colour_axes, ax=ax2, label = r'$\langle N_{\rm{sys}} \rangle$ ')
        ax2.legend()
        plt.savefig('figures/binary_hierarchical/sys_formation_N_plot.pdf', dpi=300, bbox_inches='tight')
        #p0 = (1, 1, 10**-5)
        #params, cv = scipy.optimize.curve_fit(log_fit, ini_pop[bin_formed>0], (bsys_time[bin_formed>0]), p0, maxfev = 10000, method = 'trf')
        #slope, beta, yint = params

        #print('Slope : ', pb[0])
        #print('y-int : ', pb[1])

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

    def GW_emissions(self):
        
        GW_calcs = gw_calcs()
        plot_ini = plotter_setup()

        fig = plt.figure(figsize=(8, 6))
        gs = fig.add_gridspec(2, 2,  width_ratios=(4, 2), height_ratios=(2, 4),
                              left=0.1, right=0.9, bottom=0.1, top=0.9,
                              wspace=0.05, hspace=0.05)
        ax = fig.add_subplot(gs[1, 0])
        ax1 = fig.add_subplot(gs[0, 0], sharex=ax)
        ax2 = fig.add_subplot(gs[1, 1], sharey=ax)

        hardb_idx = [[ ], [ ]]
        softb_idx = [[ ], [ ]]
        hardt_idx = [[ ], [ ]]
        softt_idx = [[ ], [ ]]
        for int_ in range(1):
            bhidx = np.where(self.hard_bin[int_] > 0)[0]
            bsidx = np.where(self.hard_bin[int_] < 0)[0]
            thidx = np.where(self.hard_ter[int_] > 0)[0]
            tsidx = np.where(self.hard_ter[int_] < 0)[0]
            hardb_idx[int_].append(bhidx)
            softb_idx[int_].append(bsidx)
            hardt_idx[int_].append(thidx)
            softt_idx[int_].append(tsidx)

        ####### CALCULATE FOR TERTIARY, HARD BINARY, SOFT BINARY AND ALL
        for int_ in range(1):
            GW_calcs.scatter_hist(self.GW_freqbin[int_], self.GW_strainbin[int_],
                                  self.GW_freqter[int_], self.GW_strainter[int_],
                                  ax, ax1, ax2, 'Binary', 'Tertiary')
            ax.set_xlabel(r'$\log_{10}f$ [Hz]')
            ax.set_ylabel(r'$\log_{10}h$')
            ax1.set_title(str(self.integrator[int_]))
            plot_ini.tickers(ax, 'plot')
            plot_ini.tickers(ax1, 'plot')
            plot_ini.tickers(ax2, 'plot')
            ax.set_ylim(-30, -12.2)
            ax.set_xlim(-12.5, 0.1)
            plt.savefig('figures/binary_hierarchical/'+str(self.integrator[int_])+'GW_freq_strain_maximise_diagram.png', dpi = 500, bbox_inches='tight')
            plt.clf()
"""

print('...sustainable_bintert_plotters...')
cst = sustainable_sys()
cst.system_formation_plotter()
cst.GW_emissions()"""