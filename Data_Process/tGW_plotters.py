from amuse.lab import *
from spatial_plotters import *
from file_logistics import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import warnings

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
        
        np.seterr(divide='ignore')
        warnings.filterwarnings("ignore", category=RuntimeWarning) 

        H0 = 73.04 #Taken from arXiv:2112.04510 in km/s/Mpc
        self.tH = (H0*(3.2408*10**-20))**-1 * 1/(3600*365*24*10**6)

        Hermite_data = glob.glob('data/Hermite/particle_trajectory/*')
        GRX_data = glob.glob('data/GRX/particle_trajectory/*')
        filename = [natsort.natsorted(Hermite_data), natsort.natsorted(GRX_data)] 
    
        self.pop = [[ ], [ ]]
        self.close_enc = [[ ], [ ]]
        self.mass_SMBH = [[ ], [ ]]
        self.mass_bin = [[ ], [ ]]

        self.semi_SMBH_avg = [[ ], [ ]]
        self.semi_SMBH_min = [[ ], [ ]]
        self.semi_bin_avg = [[ ], [ ]]
        self.semi_bin_min = [[ ], [ ]]
        self.semi_SMBH_fin = [[ ], [ ]]
        self.semi_SMBH_ini = [[ ], [ ]]
        self.ecc_SMBH_fin = [[ ], [ ]]
        self.ecc_SMBH_ini = [[ ], [ ]]
        self.time_steps = [[ ], [ ]]

        self.ecc_SMBH_avg = [[ ], [ ]]
        self.ecc_SMBH_min = [[ ], [ ]]
        self.ecc_bin_avg = [[ ], [ ]]
        self.ecc_bin_min = [[ ], [ ]]

        for int_ in range(2):
            for file_ in range(len(filename[int_])):
                with open(filename[int_][file_], 'rb') as input_file:
                    print('Reading file', file_, ': ', filename[int_][file_])
                    data = pkl.load(input_file)

                    avg_mass = 0
                    for parti_ in range(np.shape(data)[0]-1):
                        avg_mass +=  data.iloc[parti_+1][0][1].value_in(units.kg)
                    avg_mass /= (np.shape(data)[0])

                    for parti_ in range(np.shape(data)[0]): #Iterate over all particles present
                        semi_SMBH_temp = []
                        ecc_SMBH_temp = []
                        min_GW_SMBH = []
                        min_GW_bin = []
                        mass_bin = []
                        mass1 = data.iloc[parti_][0][1]
                        SMBH_sys_mass = data.iloc[0][0][1]
                        self.mass_SMBH[int_].append([mass1, SMBH_sys_mass])

                        Nclose = 0
                        bin_key = []
                        bin_sem = []
                        bin_ecc = []

                        if parti_ == 0:
                            pass
                        else:
                            self.pop[int_].append(np.shape(data)[0])
                            self.semi_SMBH_fin[int_].append(data.iloc[parti_][-2][7][0])
                            self.semi_SMBH_ini[int_].append(data.iloc[parti_][2][7][0]) 
                            self.ecc_SMBH_fin[int_].append(data.iloc[parti_][-2][8][0])
                            self.ecc_SMBH_ini[int_].append(data.iloc[parti_][2][8][0])
                            self.time_steps[int_].append(np.shape(data)[1])

                            for col_ in range(np.shape(data)[1]-1):
                                if data.iloc[parti_][col_][8][0] < 1:
                                    semi_SMBH_temp.append(data.iloc[parti_][col_][7][0].value_in(units.pc))
                                    ecc_SMBH_temp.append(data.iloc[parti_][col_][8][0])
                                    val = (data.iloc[parti_][col_][7][0].value_in(units.pc))**4 * (1-data.iloc[parti_][col_][8][0]**2)**3.5
                                    min_GW_SMBH.append(val)
                                
                                if abs(data.iloc[parti_][col_][8][1]) < 1 \
                                    and data.iloc[parti_][col_][7][1] < 0.02 | units.parsec:
                                    for part_ in range(np.shape(data)[0]):
                                        if data.iloc[part_][0][0] == data.iloc[parti_][col_][6][1]:
                                            mass2 = data.iloc[part_][0][1]
                                            bin_BE = ((constants.G*mass1*mass2)/(2*data.iloc[parti_][col_][7][1])).value_in(units.J) 
                                            if bin_BE > (150000)**2*(1+mass1/mass2):  #Hard binary conditions based on Quinlan 1996b
                                                bin_key.append(data.iloc[parti_][col_][6][1])
                                                bin_sem.append(data.iloc[parti_][col_][7][1].value_in(units.pc))
                                                bin_ecc.append(data.iloc[parti_][col_][8][1])
                                                mass_bin.append([mass1.value_in(units.MSun), mass2.value_in(units.MSun)])
                                                val = (data.iloc[parti_][col_][7][1].value_in(units.pc))**4 * (1-data.iloc[parti_][col_][8][1]**2)**3.5
                                                min_GW_bin.append(val)

                                if data.iloc[parti_][col_][-1] < 1e-2:
                                    Nclose += 1
                            self.close_enc[int_].append(Nclose)

                            mass_bin = np.asarray(mass_bin)
                            bin_in = 0
                            if len(bin_key) > 1:
                                bin_sem = np.asarray(bin_sem)
                                bin_ecc = np.asarray(bin_ecc)
                                for rep_ in range(len(bin_key)):

                                    if bin_key[rep_] != bin_key[rep_-1]:
                                        self.ecc_bin_avg[int_].append(np.mean(bin_ecc[bin_in:rep_+1]))
                                        self.semi_bin_avg[int_].append(np.mean(bin_sem[bin_in:rep_+1]))

                                        idx = np.nanargmin(min_GW_bin)
                                        self.ecc_bin_min[int_].append(bin_ecc[idx])
                                        self.semi_bin_min[int_].append(bin_sem[idx])
                                        self.mass_bin[int_].append([mass_bin[idx][0], mass_bin[idx][1]])

                                        bin_in = rep_

                                    if rep_ == len(bin_key)-1:
                                        self.ecc_bin_avg[int_].append(np.mean(bin_ecc[bin_in:rep_+1]))
                                        self.semi_bin_avg[int_].append(np.mean(bin_sem[bin_in:rep_+1]))

                                        idx = np.nanargmin(min_GW_bin)
                                        self.ecc_bin_min[int_].append(bin_ecc[idx])
                                        self.semi_bin_min[int_].append(bin_sem[idx])
                                        self.mass_bin[int_].append([mass_bin[idx][0], mass_bin[idx][1]])

                                        bin_in = rep_

                            if len(bin_key) == 1:
                                self.semi_bin_avg[int_].append(bin_sem[0])
                                self.semi_bin_min[int_].append(bin_sem[0])
                                self.ecc_bin_avg[int_].append(bin_ecc[0])
                                self.ecc_bin_min[int_].append(bin_ecc[0])
                                self.mass_bin[int_].append([mass_bin[0][0], mass_bin[0][1]])

                            semi_SMBH_temp = np.asarray(semi_SMBH_temp)
                            ecc_SMBH_temp = np.asarray(ecc_SMBH_temp)
                            idx_SMBH = np.nanargmin(min_GW_SMBH)
                            
                            self.semi_SMBH_avg[int_].append(np.mean(semi_SMBH_temp))
                            self.ecc_SMBH_avg[int_].append(np.mean(ecc_SMBH_temp))
                            self.semi_SMBH_min[int_].append(semi_SMBH_temp[idx_SMBH])
                            self.ecc_SMBH_min[int_].append(ecc_SMBH_temp[idx_SMBH])

            
            self.ecc_SMBH_avg[int_] = np.asarray(self.ecc_SMBH_avg[int_])
            self.ecc_SMBH_min[int_] = np.asarray(self.ecc_SMBH_min[int_])
            self.semi_SMBH_avg[int_] = np.asarray(self.semi_SMBH_avg[int_]) * 1 | units.parsec
            self.semi_SMBH_min[int_] = np.asarray(self.semi_SMBH_min[int_]) * 1 | units.parsec

            self.semi_bin_avg[int_] = np.asarray(self.semi_bin_avg[int_]) * 1 | units.parsec
            self.semi_bin_min[int_] = np.asarray(self.semi_bin_min[int_]) * 1 | units.parsec
            self.ecc_bin_avg[int_] = np.asarray(self.ecc_bin_avg[int_])
            self.ecc_bin_min[int_]= np.asarray(self.ecc_bin_min[int_])

            self.mass_bin[int_] = np.asarray(self.mass_bin[int_]) * 1 | units.MSun
            self.pop[int_] = np.asarray(self.pop[int_])
            self.pop[int_][self.pop[int_] %10 != 0] -= 1
            self.pop[int_] = np.asarray(self.pop[int_])

    def forecast_interferometer(self, ax, mass_arr):
        """
        Function to plot the LISA and aLIGO frequency range in Ge a vs. (1-e) plots
        """
        
        self.ecc_range = np.linspace(0.0001, 0.9999999, 50)

        self.LIGO_semimaj_min = self.peak_frequency_GW(self.ecc_range[1:], 10 | units.Hz, mass_arr[0][0], mass_arr[0][1])
        self.LIGO_semimaj = self.peak_frequency_GW(self.ecc_range[1:], 200 | units.Hz, mass_arr[0][0], mass_arr[0][1])
        self.LIGO_semimaj_max = self.peak_frequency_GW(self.ecc_range[1:], 10000 | units.Hz, mass_arr[0][0], mass_arr[0][1])

        self.LISA_semimaj_min = self.peak_frequency_GW(self.ecc_range[1:], 1e-4 | units.Hz, mass_arr[0][0], mass_arr[0][1])
        self.LISA_semimaj = self.peak_frequency_GW(self.ecc_range[1:], 1e-2 | units.Hz, mass_arr[0][0], mass_arr[0][1])
        self.LISA_semimaj_max = self.peak_frequency_GW(self.ecc_range[1:], 1 | units.Hz, mass_arr[0][0], mass_arr[0][1])

        self.ecc_range = [np.log(1-i) for i in self.ecc_range[1:]]
        self.text_angle = np.degrees(np.arctan((self.ecc_range[-4]-self.ecc_range[-3])/(self.LIGO_semimaj[-4]-self.LIGO_semimaj[-3]))) - 4

        ax.plot(self.LIGO_semimaj_min, self.ecc_range, linestyle = ':', color = 'black')
        ax.plot(self.LIGO_semimaj, self.ecc_range, linestyle = '-.', color = 'black')
        ax.plot(self.LIGO_semimaj_max, self.ecc_range, linestyle = ':', color = 'black')
        ax.fill_between(np.append(self.LIGO_semimaj_min, self.LIGO_semimaj_max[::-1]), 
                        np.append(self.ecc_range[:], self.ecc_range[::-1]), alpha = 0.2, color = 'blue')
        ax.plot(self.LISA_semimaj_min, self.ecc_range, linestyle = ':', color = 'black')
        ax.plot(self.LISA_semimaj, self.ecc_range, linestyle = '-.', color = 'black')
        ax.plot(self.LISA_semimaj_max, self.ecc_range, linestyle = ':', color = 'black')
        ax.fill_between(np.append(self.LISA_semimaj_min, self.LISA_semimaj_max[::-1]), 
                        np.append(self.ecc_range[:], self.ecc_range[::-1]), alpha = 0.2, color = 'red')

        return ax

    def gw_amplitude(self, semi, ecc, m1, m2):
        """
        Use of eqn (3) Matsubayashi et al. 2004 to calculate the approximate GW strain during the inspiral phase.
        
        Inputs:
        semi:   The semi-major axes of the system
        ecc:    The eccentricity of the binary system
        m1/m2:  The individual mass components of the binary system
        """

        dist = 8 | units.kpc  # Taken from [https://imagine.gsfc.nasa.gov/features/cosmic/milkyway_info.html]
        freq = np.sqrt(constants.G*(m1+m2)/semi**3) * (1+ecc)**1.1954 * (np.pi*(1-ecc**2)**1.5)**-1  #Frequency based on eqn 43 of Samsing et al. 2014
        hinsp = (np.sqrt(32/5) * (np.pi)**(2/3) * (constants.G)**(5/3) * constants.c**(-4) * m1*m2*(m1+m2)**(-1/3) * freq**(2/3) * dist**-1)
        hinsp = hinsp * (1 | units.kg)**5.551115123125783e-17 * (1 | units.s)**1.1102230246251565e-16

        return hinsp

    def gw_calc(self, semi, ecc, m1, m2):
        """
        Function to calculate the GW timescale
        
        Inputs:
        semi:    The semi-major axis of the binary
        ecc:     The eccentricity of the binary
        m1/m2:   The binary component masses
        outputs: The gravitational wave timescale
        """

        red_mass = (m1*m2)/(m1+m2)
        tot_mass = m1 + m2
        tgw = (5/256) * (constants.c)**5/(constants.G**3)*(semi**4*(1-ecc**2)**3.5)/(red_mass*tot_mass**2)
        return tgw

    def peak_frequency_GW(self, ecc_arr, freq_val, m1, m2):
        """
        Function to get constant frequency curves based on eqn. 43 of Samsing et al. 2014.
        Frequency values are based on Samsing et al. 2014 and correspond to LIGO (200 Hz) with range 10 < f < 10 000 [https://link.aps.org/doi/10.1103/PhysRevD.93.112004]
        LISA (1e-2 Hz) peak sensitivity with range 1e-4 < f < 1 [https://lisa.nasa.gov/]
        Inspiral Frequency Phase
        
        Inputs:
        ecc_arr:  The eccentricity array
        freq_val: The constant frequency wishing to plot
        m1/m2:    The binary mass
        output:   Constant frequency line for various semi-maj/eccentricity of binaries
        """
        
        term1 = np.sqrt(constants.G*(m1+m2))/np.pi
        semi_maj = [np.log10(((term1 * (1+i)**1.1954/(1-i**2)**1.5 * freq_val**-1)**(2/3)).value_in(units.pc)) for i in ecc_arr]
        return semi_maj

    def amp_tgw_plotter(self):
        """
        Function which plots the amplitude histogram for all sets of particle 
        using Gultekin et al. 2005 eqn 1.
        """

        plot_ini = plotter_setup()
        tgw_amp_avg_IMBH = [[ ], [ ]]
        tgw_amp_avg_SMBH = [[ ], [ ]]
        integrator = ['Hermite', 'GRX']
        for int_ in range(2):
            for i in range(len(self.semi_bin_avg[int_])):
                avg_amp = self.gw_amplitude(self.semi_bin_avg[int_][i], self.ecc_bin_avg[int_][i], 
                                            self.mass_bin[int_][i][0], self.mass_bin[int_][i][1])
                tgw_amp_avg_IMBH[int_].append(avg_amp)
            for i in range(len(self.semi_SMBH_avg[int_])):
                avg_amp = self.gw_amplitude(self.semi_SMBH_avg[int_][i], self.ecc_SMBH_avg[int_][i], 
                                            self.mass_SMBH[int_][i][0], self.mass_SMBH[int_][i][1])
                tgw_amp_avg_SMBH[int_].append(avg_amp)
            tgw_amp_avg_IMBH[int_] = np.asarray(tgw_amp_avg_IMBH[int_])
            tgw_amp_avg_SMBH[int_] = np.asarray(tgw_amp_avg_SMBH[int_])
                
            fig, ax = plt.subplots()
            ax.set_title('GW Strain Histogram')
            n, bins, patches = ax.hist(np.log10(tgw_amp_avg_IMBH[int_]), 20, histtype = 'step', density = True, color='black')
            n, bins, patches = ax.hist(np.log10(tgw_amp_avg_IMBH[int_]), 20, color='black', density = True, alpha = 0.3, label = r'IMBH-IMBH')
            n, bins, patches = ax.hist(np.log10(tgw_amp_avg_SMBH[int_]), 20, histtype = 'step', density = True, color='purple')
            n, bins, patches = ax.hist(np.log10(tgw_amp_avg_SMBH[int_]), 20, color='purple', density = True, alpha = 0.3, label = r'IMBH-SMBH')
            ax.set_yscale('log')
            ax.set_xlabel(r'$\log_{10}h$')
            ax.set_ylabel(r'Occurence')
            ax.set_title('Gravitational Strain Histogram')
            ax.text(-16, 0.6, r'$r_{merger} = 8$ kpc')
            plot_ini.tickers(ax, 'plot')
            ax.legend()
            plt.savefig('figures/gravitational_waves/GWstrain_IMBH_histogram'+str(integrator[int_])+'.pdf', dpi = 300, bbox_inches='tight')

            with open('figures/gravitational_waves/'+str(integrator[int_])+'redshift_strain.txt', 'w') as file:
                file.write('\nFor z ~ 3.5, IMBH-IMBH strain which is inversely proportional to distance becomes: ' + str(np.mean(tgw_amp_avg_IMBH[int_]) * (8000/6025) * 10**-6))
                file.write('\nFor z ~ 3.5, SMBH-IMBH strain which is inversely proportional to distance becomes: ' + str(np.mean(tgw_amp_avg_SMBH[int_]) * (8000/6025) * 10**-6))

    def SMBH_tgw_plotter(self):
        """
        Function which manipulates data to extract:
        - Minimum and average tgw w.r.t SMBH for each particle
        - Minimum and average semi-major axis w.r.t SMBH for each particle
        - Minimum and average eccentricity w.r.t SMBH for each particle
        It plots a histogram to convey the evolution of this binary over the course of the simulation.
        It also plots a vs. (1-e) in both a histogram and scatter of the minimum to compare with LISA and LIGO.
        """

        plot_init = plotter_setup()
        xmin_evol = []
        tgw_SMBH_fin = [[ ], [ ]]
        tgw_SMBH_ini = [[ ], [ ]]
        tgw_SMBH_min = [[ ], [ ]]
        tgw_SMBH_iratio = [[ ], [ ]]
        tgw_SMBH_fratio = [[ ], [ ]]
        tgw_ratio = [[ ], [ ]]
        ecc_SMBH_conc = [[ ], [ ]]
        sem_SMBH_conc = [[ ], [ ]]
        ecc_SMBH_min = [[ ], [ ]]
        ecc_SMBH_avg = [[ ], [ ]]
        semi_SMBH_min = [[ ], [ ]]
        semi_SMBH_avg = [[ ], [ ]]

        ############## PLOTTING OF TGW_F/TGW_0 ##############
        fig, ax = plt.subplots(figsize=(12,6), nrows=1, ncols=2)
        gs = gridspec.GridSpec(1, 2)
        ax1 = plt.subplot(gs[0,0])
        ax2 = plt.subplot(gs[0,1])
        ax1.set_title('Hermite')
        ax2.set_title('GRX')
        ymin_evol = [ ] ; ymax_evol = [ ]
        xmin_evol = [ ] ; xmax_evol = [ ]
        ymin_ecc_sem = [ ]
        xmin_ecc_sem = [ ] ; xmax_ecc_sem = [ ]
        for int_ in range(2):
            for i in range(len(self.semi_SMBH_avg[int_])):
                grav_f_time = self.gw_calc(self.semi_SMBH_fin[int_][i], self.ecc_SMBH_fin[int_][i], 
                                           self.mass_SMBH[int_][i][0], self.mass_SMBH[int_][i][1]).value_in(units.Myr)
                grav_i_time = self.gw_calc(self.semi_SMBH_ini[int_][i], self.ecc_SMBH_ini[int_][i], 
                                           self.mass_SMBH[int_][i][0], self.mass_SMBH[int_][i][1]).value_in(units.Myr)
                grav_min_time = self.gw_calc(self.semi_SMBH_min[int_][i], self.ecc_SMBH_min[int_][i], 
                                             self.mass_SMBH[int_][i][0], self.mass_SMBH[int_][i][1]).value_in(units.Myr)
                tgw_SMBH_fin[int_].append(grav_f_time)
                tgw_SMBH_ini[int_].append(grav_i_time)
                tgw_SMBH_min[int_].append(grav_min_time)
                
            tgw_SMBH_fin_arr = np.asarray(tgw_SMBH_fin[int_])
            tgw_SMBH_ini_arr = np.asarray(tgw_SMBH_ini[int_])
            tgw_SMBH_min_raw_arr = np.asarray(tgw_SMBH_min[int_])
            ecc_SMBH_avg_raw_arr = np.asarray(self.ecc_SMBH_avg[int_])
            ecc_SMBH_min_raw_arr = np.asarray(self.ecc_SMBH_min[int_])
            semi_SMBH_avg_raw_arr = np.asarray([i.value_in(units.parsec) for i in self.semi_SMBH_avg[int_]])
            semi_SMBH_min_raw_arr = np.asarray([i.value_in(units.parsec) for i in self.semi_SMBH_min[int_]])

            ecc_SMBH_avg[int_] = np.asarray([np.log10(1-i) for i in ecc_SMBH_avg_raw_arr[(ecc_SMBH_avg_raw_arr< 1) & (semi_SMBH_avg_raw_arr > 0)]])
            semi_SMBH_avg[int_] = np.asarray([np.log(i) for i in semi_SMBH_avg_raw_arr[(ecc_SMBH_avg_raw_arr < 1) & (semi_SMBH_avg_raw_arr > 0)]])
            ecc_SMBH_min[int_] = np.asarray([np.log10(1-i) for i in ecc_SMBH_min_raw_arr[(ecc_SMBH_min_raw_arr < 1) & (semi_SMBH_min_raw_arr > 0)]])
            semi_SMBH_min[int_] = np.asarray([np.log10(i) for i in semi_SMBH_min_raw_arr[(ecc_SMBH_min_raw_arr < 1) & (semi_SMBH_min_raw_arr > 0)]])
            tgw_SMBH_min[int_] = np.asarray(tgw_SMBH_min_raw_arr[(ecc_SMBH_min_raw_arr < 1) & (semi_SMBH_min_raw_arr > 0)])

            ecc_SMBH_conc[int_] = np.concatenate((ecc_SMBH_avg[int_], ecc_SMBH_min[int_]), axis = None)
            sem_SMBH_conc[int_] = np.concatenate((semi_SMBH_avg[int_], semi_SMBH_min[int_]), axis = None)

            tgw_SMBH_fratio[int_] = np.asarray([i / self.tH for i in tgw_SMBH_fin_arr])
            tgw_SMBH_iratio[int_] = np.asarray([i / self.tH for i in tgw_SMBH_ini_arr])
            tgw_ratio[int_] = np.asarray([i / j for i, j in zip(tgw_SMBH_fin_arr, tgw_SMBH_ini_arr)])

            xmin_evol.append(np.nanmin(np.log10(tgw_SMBH_iratio[int_])))
            xmax_evol.append(np.nanmax(np.log10(tgw_SMBH_iratio[int_])))
            ymin_evol.append(np.nanmin(np.log10(tgw_SMBH_fratio[int_])))
            ymax_evol.append(np.nanmax(np.log10(tgw_SMBH_fratio[int_])))

            xmin_ecc_sem.append(np.nanmin(semi_SMBH_min[int_]))
            xmax_ecc_sem.append(np.nanmax(semi_SMBH_min[int_]))
            ymin_ecc_sem.append(np.nanmin(ecc_SMBH_min[int_]))

        xlist = np.linspace(0, 1.1*max(xmax_evol))
        iter = -1
        for ax_ in [ax1, ax2]:
            iter += 1
            plot_init.tickers(ax_, 'hist')
            ax_.set_ylabel(r'$log_{10}(t_{\rm{GW},f}/t_{\rm{H}})$')
            ax_.set_xlabel(r'$log_{10}(t_{\rm{GW},0}/t_{\rm{H}})$')
            ax_.plot(xlist, xlist, linestyle = ':', color = 'white')
            ax_.set_xlim(0.9*min(xmin_evol), 1.1*max(xmax_evol))
            ax_.set_ylim(0.9*min(ymin_evol), 1.1*max(ymax_evol))
            bin2d_sim, xedg, xedg, image = ax_.hist2d(np.log10(tgw_SMBH_iratio[iter]), np.log10(tgw_SMBH_fratio[iter]), \
                                                              bins=(150,150), range=([0, 1.1*max(xmax_evol)], [0, 1.1*max(ymax_evol)]))
            bin2d_sim /= np.max(bin2d_sim)
            extent = [0, 1.1*max(xmax_evol), 0, 1.1*max(ymax_evol)]
            contours = ax_.imshow(np.log10(bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
        plt.savefig('figures/gravitational_waves/tgw_SMBH_evol.pdf', dpi=300, bbox_inches='tight')

        ############## PLOTTING OF a vs. (1-e) FOR SMBH ##############
        integrator = ['Hermite', 'GRX']
        for int_ in range(2):
            fig = plt.figure(figsize=(12.5, 6))
            ax1 = fig.add_subplot(121)
            ax2 = fig.add_subplot(122)
            for ax_ in [ax1, ax2]:
                ax_.set_ylabel(r'$\log_{10}(1-e)$')
                ax_.set_xlabel(r'$\log_{10} a$ [pc]')
                ax_.set_ylim(1.1*min(ymin_ecc_sem), 0)
                ax_.set_xlim(-8, 1.1*max(xmax_ecc_sem))
            bin2d_sim, xedg, xedg, image = ax2.hist2d(sem_SMBH_conc[int_], ecc_SMBH_conc[int_], bins=(100,100), 
                                                      range=([1.1*min(xmin_ecc_sem), 1.1*max(xmax_ecc_sem)], [1.1*min(ymin_ecc_sem), 0]))
            idx_merge = np.where(tgw_SMBH_min[int_] <= (self.tH/10**6))
            idx_stab = np.where(tgw_SMBH_min[int_] > (self.tH/10**6))
            
            norm_min = (min(np.nanmin(tgw_SMBH_min[0]), np.nanmin(tgw_SMBH_min[0]),
                            np.nanmin(tgw_SMBH_min[1]), np.nanmin(tgw_SMBH_min[1])))
            norm_max = (max(np.nanmax(tgw_SMBH_min[0]), np.nanmax(tgw_SMBH_min[0]),
                            np.nanmax(tgw_SMBH_min[1]), np.nanmax(tgw_SMBH_min[1])))
            normalise = plt.Normalize((norm_min), (norm_max))

            bin2d_sim /= np.max(bin2d_sim)
            extent = [1.1*min(xmin_ecc_sem), 1.1*max(xmax_ecc_sem), 1.1*min(ymin_ecc_sem), 0]
            contours = ax2.imshow(bin2d_sim, extent = extent, aspect='auto', origin = 'upper')
            plot_init.tickers(ax2, 'hist')
            plot_init.tickers(ax1, 'plot')
            colour_axes = ax1.scatter(semi_SMBH_min[int_][idx_merge], ecc_SMBH_min[int_][idx_merge], norm = normalise,
                                      c = np.log10(tgw_SMBH_min[int_][idx_merge]), marker = 'X', edgecolors='black')
            colour_axes = ax1.scatter(semi_SMBH_min[int_][idx_stab], ecc_SMBH_min[int_][idx_stab], norm = normalise,
                                      c = np.log10(tgw_SMBH_min[int_][idx_stab]), edgecolors='black')
            self.forecast_interferometer(ax1, self.mass_SMBH[int_])
            plt.colorbar(colour_axes, ax = ax1, label = r'$\log_{10} \langle t_{GW}\rangle$ [Myr]')
            ax1.text(-6.87, -6, r'$f = 200$ Hz', verticalalignment = 'center', fontsize ='xx-small', rotation=self.text_angle-2, color = 'black')
            ax1.text(-4., -6, r'$f = 10^{-2}$ Hz', verticalalignment = 'center', fontsize ='xx-small', rotation=self.text_angle-2, color = 'black')
            plt.savefig('figures/gravitational_waves/ecc_semi_SMBH'+str(integrator[int_])+'_histogram.pdf', dpi=300, bbox_inches='tight')

        ############## PLOTTING OF N vs. <tGW> FOR SMBH ##############

        fig = plt.figure(figsize=(14, 6))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax1.set_title('Hermite')
        ax2.set_title('GRX')
        ax_ = [ax1, ax2]
        ymin = [ ] ; ymax = [ ]
        norm_min = min(np.nanmin(self.close_enc[0]), np.nanmin(self.close_enc[1]))
        norm_max = max(np.nanmax(self.close_enc[1]), np.nanmax(self.close_enc[1]))
        normalise = plt.Normalize(norm_min, norm_max)
        for int_ in range(2):   
            ini_pop = np.unique(self.pop[int_])
            avg_tgw = np.empty(len(ini_pop))
            close_enc = np.empty(len(ini_pop))
            iter = -1   
            for pop_ in ini_pop:
                iter +=1
                idx = np.where(self.pop[int_] == pop_)[0]
                avg_tgw[iter] = np.mean(tgw_ratio[int_][idx])
                close_enc[iter] = np.mean(np.asarray(self.close_enc[int_])[idx])
            colour_axes = ax_[int_].scatter(ini_pop, np.log10(avg_tgw), edgecolors = 'black', c = close_enc, norm = normalise)
            ax_[int_].set_xlabel(r'IMBH Population [$N$]')
            ax_[int_].set_ylabel(r'$\log_{10}\langle t_{\rm{GW},f}/t_{\rm{GW},0} \rangle$')
            plot_init.tickers_pop(ax_[int_], ini_pop)
            ymin.append(np.nanmin(np.log10(avg_tgw)))
            ymax.append(np.nanmax(np.log10(avg_tgw)))
            
        for axis in ax_:
            axis.set_ylim(1.05*min(ymin), 1.05*max(ymax))
        plt.colorbar(colour_axes, ax = ax2, label = r'# Close Encounters')
        plt.savefig('figures/gravitational_waves/pop_tgwevol_enc_plot.pdf', dpi=300, bbox_inches='tight')

        fig = plt.figure(figsize=(14, 6))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax1.set_title('Hermite')
        ax2.set_title('GRX')

        xmax = max(np.nanmax(np.log10(self.close_enc[0])), np.nanmax(np.log10(self.close_enc[1])))
        ymin = min(np.nanmin(np.log10(tgw_ratio[0])), np.nanmin(np.log10(tgw_ratio[0])))
        ymax = max(np.nanmax(np.log10(tgw_ratio[0])), np.nanmax(np.log10(tgw_ratio[0])))
        iter = -1
        for ax_ in [ax1, ax2]:
            iter += 1
            bin2d_sim, xedg, xedg, image = ax_.hist2d(np.log10(self.close_enc[iter]), np.log10(tgw_ratio[iter]), 
                                                      bins=(50,50), range=([0, 1.1*xmax], [0.9*ymin, 1.1*ymax]))
            bin2d_sim /= np.max(bin2d_sim)
            extent = [0, 1.1*xmax, 0.9*ymin, 1.1*ymax]
            contours = ax_.imshow((bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
            ax_.set_xlabel(r'$\log_{10} N_{\rm{enc}}$')
            ax_.set_ylabel(r'$\log_{10}(t_{\rm{GW},f}/t_{\rm{GW},0})$')
            ax_.set_xlim(0, 1.1*xmax)
            plot_init.tickers(ax_, 'hist')
        plt.savefig('figures/gravitational_waves/Nenc_tgw_SMBH_evol_histogram.pdf', dpi=300, bbox_inches='tight')
        plt.clf()

        for int_ in range(2):
            rISCO = 6*constants.G*self.mass_SMBH[int_][0][1] * constants.c**-2
            freqGW = 2*np.sqrt(constants.G*self.mass_SMBH[int_][0][1]/(4*np.pi**2*rISCO**3))
            redshift2 = freqGW/(10**-4 | units.Hz) - 1

            avg_SMBH_freq = np.sqrt(constants.G*(self.mass_SMBH[int_][0][0]+self.mass_SMBH[int_][0][1])/np.mean(self.semi_SMBH_avg[int_][self.semi_SMBH_avg[int_] > 0 | units.m])**3) \
                                    * (1+np.mean(self.ecc_SMBH_avg[int_][self.ecc_SMBH_avg[int_]<1]))**1.1954 * (np.pi*(1-np.mean(self.ecc_SMBH_avg[int_][self.ecc_SMBH_avg[int_] < 1])**2)**1.5)**-1
            avg_strain = self.gw_amplitude(np.mean(self.semi_SMBH_avg[int_]), np.mean(self.ecc_SMBH_avg[int_]), self.mass_SMBH[int_][0][0], self.mass_SMBH[int_][0][1])
            redshift_avg = avg_SMBH_freq / (10**-4 | units.Hz)  - 1

            avg_SMBH_freq = np.sqrt(constants.G*(self.mass_SMBH[int_][0][0]+self.mass_SMBH[int_][0][1])/np.mean(self.semi_SMBH_min[int_][self.semi_SMBH_min[int_] >0 | units.m])**3) \
                                    * (1+np.mean(self.ecc_SMBH_min[int_][self.ecc_SMBH_min[int_] < 1]))**1.1954 * (np.pi*(1-np.mean(self.ecc_SMBH_min[int_][self.ecc_SMBH_min[int_] < 1])**2)**1.5)**-1
            redshift_min = avg_SMBH_freq / (10**-4 | units.Hz)  - 1

            with open('figures/gravitational_waves/'+str(integrator[int_])+'redshift_frequency.txt', 'w') as file:
                file.write('\nFrequency outcomes for '+str(integrator[int_]))
                file.write('\nAt r ~ rISCO the frequency takes a different expression.')
                file.write('\nTypical frequency becomes: '+str(freqGW.in_(units.Hz)))
                file.write('\nMaking them observable up to z < '+str(redshift2))
                file.write('\nGiven the avg. SMBH frequency '+str(avg_SMBH_freq.in_(units.Hz))+' and the minimum LISA frequency 1e-4 Hz.')
                file.write('\nThen to detect SMBH - IMBH merger with a supposed avg. GW amplitude '+str(avg_strain))
                file.write('\nthe cluster needs to be at a distance z < '+str(redshift_avg))
                file.write('\nor in the best case scenario z < '+str(redshift_min))
                file.write('\n========================================')

    def bin_tgw_plotter(self):
        """
        Function which detects all binaries and plots their eccentricity vs. semi-major axis
        """

        plot_init = plotter_setup()

        tgw_bin_avg = [[ ], [ ]]
        tgw_bin_min = [[ ], [ ]]
        semi_bin_min = [[ ], [ ]]
        semi_bin_avg = [[ ], [ ]]
        semi_bin_conc = [[ ], [ ]]
        ecc_bin_min = [[ ], [ ]]
        ecc_bin_avg = [[ ], [ ]]
        ecc_bin_conc = [[ ], [ ]]
        BE = [[ ], [ ]]
        KE = [[ ], [ ]]
        idx_hard = [[ ], [ ]]
        idx_soft = [[ ], [ ]]
        for int_ in range(2):
            for i in range(len(self.semi_bin_avg[int_])):
                grav_bin_avg_timescale = self.gw_calc(self.semi_bin_avg[int_][i], self.ecc_bin_avg[int_][i], 
                                                      self.mass_bin[int_][i][0], self.mass_bin[int_][i][1]).value_in(units.Myr)
                grav_bin_min_timescale = self.gw_calc(self.semi_bin_min[int_][i], self.ecc_bin_min[int_][i], 
                                                      self.mass_bin[int_][i][0], self.mass_bin[int_][i][1]).value_in(units.Myr)
                BE[int_].append((constants.G*self.mass_bin[int_][i][0]*self.mass_bin[int_][i][0]/(2*self.semi_bin_min[int_][i])).value_in(units.J))
                tgw_bin_avg[int_].append(grav_bin_avg_timescale)
                tgw_bin_min[int_].append(grav_bin_min_timescale)

            semi_bin_avg_raw = np.asarray([i.value_in(units.parsec) for i in self.semi_bin_avg[int_]])
            semi_bin_min_raw = np.asarray([i.value_in(units.parsec) for i in self.semi_bin_min[int_]])
            ecc_bin_avg_raw = np.asarray(self.ecc_bin_avg[int_])
            ecc_bin_min_raw = np.asarray(self.ecc_bin_min[int_])
            tgw_bin_min_raw = np.asarray(tgw_bin_min[int_])
            BE[int_] = np.asarray(BE[int_])
            KE[int_] = 0.5*self.mass_bin[int_][0][0].value_in(units.kg)*(150000)**2

            ecc_bin_avg[int_] = np.asarray([np.log10(1-i) for i in ecc_bin_avg_raw[(ecc_bin_avg_raw < 1) & (semi_bin_avg_raw > 0)]])
            semi_bin_avg[int_] = np.asarray([np.log10(i) for i in semi_bin_avg_raw[(ecc_bin_avg_raw < 1) & (semi_bin_avg_raw > 0)]])
            ecc_bin_min[int_] = np.asarray([np.log10(1-i) for i in ecc_bin_min_raw[(ecc_bin_min_raw < 1) & (semi_bin_min_raw > 0)]])
            semi_bin_min[int_] = np.asarray([np.log10(i) for i in semi_bin_min_raw[(ecc_bin_min_raw < 1) & (semi_bin_min_raw > 0)]])
            tgw_bin_min[int_] = np.asarray(tgw_bin_min_raw[(ecc_bin_min_raw < 1) & (semi_bin_min_raw > 0)])
            BE[int_] = np.asarray(BE[int_][(ecc_bin_min_raw < 1) & (semi_bin_min_raw > 0)])
            ecc_bin_min[int_] = np.asarray(ecc_bin_min[int_][semi_bin_min[int_] > -5])
            tgw_bin_min[int_] = np.asarray(tgw_bin_min[int_][semi_bin_min[int_] > -5])
            BE[int_] = np.asarray(BE[int_][semi_bin_min[int_] > -5])
            semi_bin_min[int_] = np.asarray(semi_bin_min[int_][semi_bin_min[int_] > -5])
            idx_hard[int_] = np.where(BE[int_] > KE[int_])[0]
            idx_soft[int_] = np.where(BE[int_] < KE[int_])[0]

            ecc_bin_conc[int_] = np.concatenate((ecc_bin_avg[int_], ecc_bin_min[int_]), axis = None)
            semi_bin_conc[int_] = np.concatenate((semi_bin_avg[int_], semi_bin_min[int_]), axis = None)        
        
        ############### PLOTTING OF a vs. (1-e) FOR BIN ##############
        #################### DO THE SAME FOR GRX #####################

        integrator = ['Hermite', 'GRX']
        
        xmin = min(np.nanmin(semi_bin_conc[0]), np.nanmin(semi_bin_conc[1]))
        xmax = max(np.nanmax(semi_bin_conc[0]), np.nanmax(semi_bin_conc[1]))
        ymin = min(np.nanmin(ecc_bin_conc[0]), np.nanmin(ecc_bin_conc[1]))
        cmin = min(np.nanmin(tgw_bin_min[0]), np.nanmin(tgw_bin_min[1]))
        cmax = max(np.nanmax(tgw_bin_min[0]), np.nanmax(tgw_bin_min[1]))
        norm_min = np.log10(cmin)
        norm_max = np.log10(cmax)
        normalise = plt.Normalize(norm_min, norm_max)

        for int_ in range(2):
            fig = plt.figure(figsize=(15, 6))
            ax1 = fig.add_subplot(121)
            ax2 = fig.add_subplot(122)
            bin2d_sim, xedg, xedg, image = ax2.hist2d(semi_bin_conc[int_], ecc_bin_conc[int_], bins=(75,75), 
                                                      range=([1.1*xmin, 1.1*xmax], [1.1*ymin, 0]))
            bin2d_sim /= np.max(bin2d_sim)
            extent = [1.1*xmin, 1.1*xmax, 1.1*ymin, 0]
            contours = ax2.imshow(np.log10(bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
            for ax_ in [ax1, ax2]:
                ax_.set_ylabel(r'$\log_{10}(1-e)$')
                ax_.set_xlabel(r'$\log_{10} a$ [pc]')
                ax_.set_ylim(1.1*ymin, 0)
                ax_.set_xlim(-10, 1.1*xmax)
                
            plot_init.tickers(ax2, 'hist')
            plot_init.tickers(ax1, 'plot')
            colour_axes = ax1.scatter((semi_bin_min[int_][idx_hard[int_]]), ecc_bin_min[int_][idx_hard[int_]], 
                                       c = np.log10(tgw_bin_min[int_][idx_hard[int_]]), 
                                       norm = normalise, edgecolors='black',  marker = 'X')
            colour_axes = ax1.scatter((semi_bin_min[int_][idx_soft[int_]]), ecc_bin_min[int_][idx_soft[int_]], 
                                       c = np.log10(tgw_bin_min[int_][idx_soft[int_]]), 
                                       norm = normalise, edgecolors='black')
            self.forecast_interferometer(ax1, self.mass_bin[int_])
            plt.colorbar(colour_axes, ax = ax1, label = r'$\log_{10} \langle t_{GW}\rangle$ [Myr]')
            ax1.text(-8.5, -5, r'$f = 200$ Hz', verticalalignment = 'center', fontsize ='xx-small', rotation=self.text_angle-3, color = 'black')
            ax1.text(-5.6, -5, r'$f = 10^{-2}$ Hz', verticalalignment = 'center', fontsize ='xx-small', rotation=self.text_angle-3, color = 'black')
            plt.savefig('figures/gravitational_waves/ecc_semi_bin_histogram'+str(integrator[int_])+'.pdf', dpi=300, bbox_inches='tight')
