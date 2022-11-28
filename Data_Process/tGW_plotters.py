from amuse.lab import *
from file_logistics import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import warnings

class coupled_systems(object):
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
        print('!!!!!! WARNING THIS WILL TAKE A WHILE !!!!!!!')

        H0 = 73.04 #Taken from arXiv:2112.04510 in km/s/Mpc
        self.tH = (H0*(3.2408*10**-20))**-1 * 1/(3600*365*24*10**6)

        Hermite_data = glob.glob('data/Hermite/particle_trajectory/*')
        GRX_data = glob.glob('data/GRX/particle_trajectory/*')
        filename = [natsort.natsorted(Hermite_data), natsort.natsorted(GRX_data)] 
    
        self.pop = [[ ], [ ]]
        self.close_enc = [[ ], [ ]]
        self.mass_SMBH = [[ ], [ ]]
        self.mass_IMBH = [[ ], [ ]]

        self.semi_SMBH_avg = [[ ], [ ]]
        self.semi_SMBH_min = [[ ], [ ]]
        self.semi_IMBH_avg = [[ ], [ ]]
        self.semi_IMBH_min = [[ ], [ ]]
        self.ecc_SMBH_avg = [[ ], [ ]]
        self.ecc_SMBH_min = [[ ], [ ]]
        self.ecc_IMBH_avg = [[ ], [ ]]
        self.ecc_IMBH_min = [[ ], [ ]]
        self.time_IMBH_avg = [[ ], [ ]]
        self.time_IMBH_min = [[ ], [ ]]

        self.freq_flyby_nn = [[ ], [ ]]
        self.freq_flyby_t = [[ ], [ ]]
        self.strain_flyby_nn = [[ ], [ ]]
        self.strain_flyby_t = [[ ], [ ]]
        self.time_flyby_nn = [[ ], [ ]]
        self.time_flyby_t = [[ ], [ ]]

        self.semi_SMBH_fin = [[ ], [ ]]
        self.semi_SMBH_ini = [[ ], [ ]]
        self.ecc_SMBH_fin = [[ ], [ ]]
        self.ecc_SMBH_ini = [[ ], [ ]]

        self.tot_sim_time = [[ ], [ ]]
        self.tot_events = [[ ], [ ]]
        self.fb_nn_events = [[ ], [ ]]
        self.fb_t_events = [[ ], [ ]]

        for int_ in range(2):
            for file_ in range(len(filename[int_])):
                with open(filename[int_][file_], 'rb') as input_file:
                    print('Reading file', file_, ': ', filename[int_][file_])
                    data = pkl.load(input_file)
                    Nclose = 0
                    Nfbnn = 0
                    Nfbt = 0
                    sim_time = np.shape(data)[1]-1

                    for parti_ in range(np.shape(data)[0]):
                        print(parti_)
                        semi_SMBH = []
                        ecc_SMBH = []
                        min_GW_SMBH = []

                        IMBH_key = []
                        IMBH_sem = []
                        IMBH_ecc = []
                        min_GW_bin = []
                        mass_IMBH = []
                        bin_time = []

                        mass1 = data.iloc[parti_][0][1]
                        SMBH_sys_mass = data.iloc[0][0][1]
                        self.mass_SMBH[int_].append([mass1, SMBH_sys_mass])

                        if parti_ != 0:
                            self.pop[int_].append(10*round(0.1*np.shape(data)[0]))
                            self.semi_SMBH_fin[int_].append(data.iloc[parti_][-2][7][0])
                            self.semi_SMBH_ini[int_].append(data.iloc[parti_][2][7][0]) 
                            self.ecc_SMBH_fin[int_].append(data.iloc[parti_][-2][8][0])
                            self.ecc_SMBH_ini[int_].append(data.iloc[parti_][2][8][0])

                            for col_ in range(np.shape(data)[1]-1):
                                if data.iloc[parti_][col_][8][0] < 1:
                                    semi_SMBH.append(data.iloc[parti_][col_][7][0].value_in(units.pc))
                                    ecc_SMBH.append(data.iloc[parti_][col_][8][0])
                                    val = (data.iloc[parti_][col_][7][0].value_in(units.pc))**4 * (1-data.iloc[parti_][col_][8][0]**2)**3.5
                                    min_GW_SMBH.append(val)
                                
                                semi_major_nn = data.iloc[parti_][col_][7][1]
                                semi_major_t = data.iloc[parti_][col_][7][2]
                                ecc_nn = (data.iloc[parti_][col_][8][1])
                                ecc_t = (data.iloc[parti_][col_][8][2])

                                for part_ in range(np.shape(data)[0]):
                                    if data.iloc[part_][0][0] == data.iloc[parti_][col_][6][1]:
                                        mass2 = data.iloc[part_][0][1]
                                
                                if np.shape(data)[0] == 10: 
                                    strain_nn, freqGW_nn = self.gw_amp_freq(semi_major_nn, ecc_nn, mass1, mass2)
                                    linex = data.iloc[parti_][col_][2][0] - data.iloc[0][col_][2][0]
                                    liney = data.iloc[parti_][col_][2][1] - data.iloc[0][col_][2][1]
                                    linez = data.iloc[parti_][col_][2][2] - data.iloc[0][col_][2][2]
                                    dist_SMBH = (linex**2+liney**2+linez**2).sqrt()
                                    dist_NN = data.iloc[parti_][col_][-1]

                                    if freqGW_nn > 10**-6:
                                        self.freq_flyby_nn[int_].append(freqGW_nn)
                                        self.strain_flyby_nn[int_].append(strain_nn)
                                        self.time_flyby_nn[int_].append(col_ * 10**-3)
                                        if dist_SMBH.value_in(units.pc) == dist_NN:
                                            print('SMBH-IMBH event')
                                            Nfbnn += 1
                                        else:
                                            Nfbnn += 0.5

                                    strain_t, freqGW_t = self.gw_amp_freq(semi_major_t, ecc_t, mass1, mass2)
                                    if freqGW_t > 10**-6:
                                        self.freq_flyby_t[int_].append(freqGW_t)
                                        self.strain_flyby_t[int_].append(strain_t)
                                        self.time_flyby_t[int_].append(col_ * 10**-3)

                                if ecc_nn < 1 and semi_major_nn < 0.02 | units.parsec:
                                    IMBH_BE = ((constants.G*mass1*mass2)/(2*data.iloc[parti_][col_][7][1])).value_in(units.J) 
                                    if IMBH_BE > (150000)**2*(1+mass1/mass2):  #Hard binary conditions based on Quinlan 1996b
                                        IMBH_key.append(data.iloc[parti_][col_][6][1])
                                        IMBH_sem.append(data.iloc[parti_][col_][7][1].value_in(units.pc))
                                        IMBH_ecc.append(data.iloc[parti_][col_][8][1])
                                        mass_IMBH.append([mass1.value_in(units.MSun), mass2.value_in(units.MSun)])
                                        val = (data.iloc[parti_][col_][7][1].value_in(units.pc))**4 * (1-data.iloc[parti_][col_][8][1]**2)**3.5
                                        min_GW_bin.append(val)
                                        bin_time.append(1000*col_)

                                if data.iloc[parti_][col_][-1] < 1e-2:
                                    Nclose += 1
                            self.close_enc[int_].append(Nclose)

                            mass_IMBH = np.asarray(mass_IMBH)
                            IMBH_in = 0
                            if len(IMBH_key) > 1:
                                IMBH_sem = np.asarray(IMBH_sem)
                                IMBH_ecc = np.asarray(IMBH_ecc)
                                for rep_ in range(len(IMBH_key)):
                                    if IMBH_key[rep_] != IMBH_key[rep_-1]:
                                        self.semi_IMBH_avg[int_].append(np.mean(IMBH_sem[IMBH_in:rep_+1]))
                                        self.ecc_IMBH_avg[int_].append(np.mean(IMBH_ecc[IMBH_in:rep_+1]))
                                        self.time_IMBH_avg[int_].append(np.mean(bin_time[IMBH_in:rep_+1]))

                                        idx = np.nanargmin(min_GW_bin)
                                        self.semi_IMBH_min[int_].append(IMBH_sem[idx])
                                        self.ecc_IMBH_min[int_].append(IMBH_ecc[idx])
                                        self.mass_IMBH[int_].append([mass_IMBH[idx][0], mass_IMBH[idx][1]])
                                        self.time_IMBH_min[int_].append(bin_time[idx])

                                        IMBH_in = rep_

                                    if rep_ == len(IMBH_key)-1: #Otherwise blocks out at the last time step
                                        self.ecc_IMBH_avg[int_].append(np.mean(IMBH_ecc[IMBH_in:rep_+1]))
                                        self.semi_IMBH_avg[int_].append(np.mean(IMBH_sem[IMBH_in:rep_+1]))
                                        self.time_IMBH_avg[int_].append(np.mean(bin_time[IMBH_in:rep_+1]))

                                        idx = np.nanargmin(min_GW_bin)
                                        self.ecc_IMBH_min[int_].append(IMBH_ecc[idx])
                                        self.semi_IMBH_min[int_].append(IMBH_sem[idx])
                                        self.mass_IMBH[int_].append([mass_IMBH[idx][0], mass_IMBH[idx][1]])
                                        self.time_IMBH_min[int_].append(bin_time[idx])

                                        IMBH_in = rep_

                            if len(IMBH_key) == 1:
                                self.semi_IMBH_avg[int_].append(IMBH_sem[0])
                                self.semi_IMBH_min[int_].append(IMBH_sem[0])
                                self.ecc_IMBH_avg[int_].append(IMBH_ecc[0])
                                self.ecc_IMBH_min[int_].append(IMBH_ecc[0])
                                self.mass_IMBH[int_].append([mass_IMBH[0][0], mass_IMBH[0][1]])
                                self.time_IMBH_avg[int_].append(bin_time[0])
                                self.time_IMBH_min[int_].append(bin_time[0])

                            semi_SMBH = np.asarray(semi_SMBH)
                            ecc_SMBH = np.asarray(ecc_SMBH)
                            idx_SMBH = np.nanargmin(min_GW_SMBH)
                            
                            self.semi_SMBH_avg[int_].append(np.mean(semi_SMBH))
                            self.ecc_SMBH_avg[int_].append(np.mean(ecc_SMBH))
                            self.semi_SMBH_min[int_].append(semi_SMBH[idx_SMBH])
                            self.ecc_SMBH_min[int_].append(ecc_SMBH[idx_SMBH])

                self.tot_sim_time[int_].append(sim_time*10**3)
                self.tot_events[int_].append(Nfbnn+Nfbt)
                self.fb_nn_events[int_].append(Nfbnn)
                self.fb_t_events[int_].append(Nfbt)

            self.ecc_SMBH_avg[int_] = np.asarray(self.ecc_SMBH_avg[int_])
            self.ecc_SMBH_min[int_] = np.asarray(self.ecc_SMBH_min[int_])
            self.semi_SMBH_avg[int_] = np.asarray(self.semi_SMBH_avg[int_]) * 1 | units.parsec
            self.semi_SMBH_min[int_] = np.asarray(self.semi_SMBH_min[int_]) * 1 | units.parsec

            self.semi_IMBH_avg[int_] = np.asarray(self.semi_IMBH_avg[int_]) * 1 | units.parsec
            self.semi_IMBH_min[int_] = np.asarray(self.semi_IMBH_min[int_]) * 1 | units.parsec
            self.ecc_IMBH_avg[int_] = np.asarray(self.ecc_IMBH_avg[int_])
            self.ecc_IMBH_min[int_]= np.asarray(self.ecc_IMBH_min[int_])

            self.mass_IMBH[int_] = np.asarray(self.mass_IMBH[int_]) * 1 | units.MSun
            self.pop[int_] = np.asarray(self.pop[int_])

            self.tot_events[int_] = np.asarray(self.tot_events[int_])
            self.tot_sim_time[int_] = np.asarray(self.tot_sim_time[int_])
            self.fb_nn_events[int_] = np.asarray(self.fb_nn_events[int_])
            self.fb_t_events[int_] = np.asarray(self.fb_t_events[int_])

        print(self.semi_SMBH_avg, self.ecc_SMBH_avg)
        print(self.semi_SMBH_min, self.ecc_SMBH_min)
        print(self.semi_IMBH_avg, self.ecc_IMBH_avg)
        print(self.semi_IMBH_min, self.ecc_IMBH_min)

        print('events')
        print(self.time_flyby_nn, self.time_flyby_t, self.freq_flyby_nn, self.freq_flyby_t)
        

        integrator = ['Hermite', 'GRX']
        with open('figures/gravitational_waves/output/event_rate.txt', 'w') as file:
            for int_ in range(2):
                file.write('\nData for '+str(integrator[int_]))
                file.write('\nAverage total event rate per Myr            ' + str(np.mean(self.tot_events[int_]/self.tot_sim_time[int_] * 10**6)))
                file.write('\nAverage nearest neigh. event rate per Myr   ' + str(np.mean(self.fb_nn_events[int_]/self.tot_sim_time[int_] * 10**6)))
                file.write('\nAverage tertiary event rate per Myr         ' + str(np.mean(self.fb_t_events[int_]/self.tot_sim_time[int_] * 10**6)))
                file.write('\n========================================================================')

    def coll_radius(self, mass_arr):
        return 3 * (2*constants.G*mass_arr)/(constants.c**2)

    def forecast_interferometer(self, ax, mass_arr):
        """
        Function to plot the LISA and aLIGO frequency range in Ge a vs. (1-e) plots
        """
        
        self.ecc_range = np.linspace(0.0001, 0.9999999, 50)

        self.LIGO_semimaj_max = self.GW_freq(self.ecc_range[1:], 10000 | units.Hz, mass_arr[0][0], mass_arr[0][1])
        self.LIGO_semimaj_min = self.GW_freq(self.ecc_range[1:], 10 | units.Hz, mass_arr[0][0], mass_arr[0][1])
        self.LIGO_semimaj = self.GW_freq(self.ecc_range[1:], 200 | units.Hz, mass_arr[0][0], mass_arr[0][1])

        self.LISA_semimaj_max = self.GW_freq(self.ecc_range[1:], 1 | units.Hz, mass_arr[0][0], mass_arr[0][1])
        self.LISA_semimaj_min = self.GW_freq(self.ecc_range[1:], 1e-4 | units.Hz, mass_arr[0][0], mass_arr[0][1])
        self.LISA_semimaj = self.GW_freq(self.ecc_range[1:], 1e-2 | units.Hz, mass_arr[0][0], mass_arr[0][1])

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
    
    def gw_amp_freq(self, semi, ecc, m1, m2):
        """
        Use of eqn (3) Matsubayashi et al. 2004 to calculate the approximate GW strain.
        Frequency equation is based on Samsing et al. 2014 eqn (43). 
        
        Inputs:
        semi:   The semi-major axes of the system
        ecc:    The eccentricity of the binary system
        m1/m2:  The individual mass components of the binary system
        """

        dist = 1 | units.Mpc  # Taken from [https://imagine.gsfc.nasa.gov/features/cosmic/milkyway_info.html]
        freq = np.sqrt(constants.G*(m1+m2)/semi**3) * (1+ecc)**1.1954 * (np.pi*(1-ecc**2)**1.5)**-1  #Frequency based on eqn 43 of Samsing et al. 2014
        strain = (np.sqrt(32/5) * (np.pi)**(2/3) * (constants.G)**(5/3) * constants.c**(-4) * m1*m2*(m1+m2)**(-1/3) * freq**(2/3) * dist**-1)
        strain = strain * (1 | units.kg)**5.551115123125783e-17 * (1 | units.s)**1.1102230246251565e-16

        return strain, freq.value_in(units.Hz)

    def gw_calc(self, semi, ecc, m1, m2):
        """
        Function to calculate the GW timescale based on Peters (1964).
        
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

    def GW_freq(self, ecc_arr, freq_val, m1, m2):
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

    def strain_freq_plotter(self):
        """
        Function which plots the amplitude histogram for all sets of particle 
        using Gultekin et al. 2005 eqn 1.
        """

        plot_ini = plotter_setup()
        integrator = ['Hermite', 'GRX']
        tgw_amp_avg_IMBH = [[ ], [ ]]
        tgw_frq_avg_IMBH = [[ ], [ ]]
        tgw_amp_avg_SMBH = [[ ], [ ]]
        tgw_frq_avg_SMBH = [[ ], [ ]]

        tgw_amp_min_IMBH = [[ ], [ ]]
        tgw_frq_min_IMBH = [[ ], [ ]]
        tgw_amp_min_SMBH = [[ ], [ ]]
        tgw_frq_min_SMBH = [[ ], [ ]]
        
        fig = plt.figure(figsize=(12.5, 6))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax_ = [ax1, ax2]
        for int_ in range(2):
            for i in range(len(self.semi_IMBH_avg[int_])):
                ecc_IMBH = self.ecc_IMBH_avg[int_][i]
                sem_IMBH = self.semi_IMBH_avg[int_][i]
                avg_amp, avg_freq = self.gw_amp_freq(sem_IMBH, ecc_IMBH, self.mass_IMBH[int_][i][0], self.mass_IMBH[int_][i][1])
                tgw_amp_avg_IMBH[int_].append(avg_amp)
                tgw_frq_avg_IMBH[int_].append(avg_freq)

                ecc_IMBHm = self.ecc_IMBH_min[int_][i]
                sem_IMBHm = self.semi_IMBH_min[int_][i]
                min_amp, min_freq = self.gw_amp_freq(sem_IMBHm, ecc_IMBHm, self.mass_IMBH[int_][i][0], self.mass_IMBH[int_][i][1])
                tgw_amp_min_IMBH[int_].append(min_amp)
                tgw_frq_min_IMBH[int_].append(min_freq)

            for k in range(len(self.semi_SMBH_min[int_])):
                ecc_SMBH = self.ecc_SMBH_avg[int_][k]
                sem_SMBH = self.semi_SMBH_avg[int_][k]
                avg_amp, avg_freq = self.gw_amp_freq(sem_SMBH, ecc_SMBH, self.mass_SMBH[int_][k][0], self.mass_SMBH[int_][k][1])
                tgw_amp_avg_SMBH[int_].append(avg_amp)
                tgw_frq_avg_SMBH[int_].append(avg_freq)

                ecc_SMBHm = self.ecc_SMBH_min[int_][k]
                sem_SMBHm = self.semi_SMBH_min[int_][k]
                min_amp, min_freq = self.gw_amp_freq(sem_SMBHm, ecc_SMBHm, self.mass_SMBH[int_][k][0], self.mass_SMBH[int_][k][1])
                tgw_amp_min_SMBH[int_].append(min_amp)
                tgw_frq_min_SMBH[int_].append(min_freq)

            tgw_amp_avg_IMBH[int_] = np.asarray(tgw_amp_avg_IMBH[int_])
            tgw_frq_avg_IMBH[int_] = np.asarray(tgw_frq_avg_IMBH[int_])
            tgw_amp_avg_SMBH[int_] = np.asarray(tgw_amp_avg_SMBH[int_])
            tgw_frq_avg_SMBH[int_] = np.asarray(tgw_frq_avg_SMBH[int_])
                
            n, bins, patches = ax_[int_].hist(np.log10(tgw_amp_avg_SMBH[int_]), 20, histtype = 'step', density = True, color='black')
            n, bins, patches = ax_[int_].hist(np.log10(tgw_amp_avg_SMBH[int_]), 20, color='black', density = True, alpha = 0.3, label = r'SMBH-IMBH')
            n, bins, patches = ax_[int_].hist(np.log10(tgw_amp_avg_IMBH[int_]), 20, histtype = 'step', density = True, color='purple')
            n, bins, patches = ax_[int_].hist(np.log10(tgw_amp_avg_IMBH[int_]), 20, color='purple', density = True, alpha = 0.3, label = r'IMBH-IMBH')
            ax_[int_].set_yscale('log')
            ax_[int_].set_xlabel(r'$\log_{10}h$')
            ax_[int_].set_title(str(integrator[int_]))
            plot_ini.tickers(ax_[int_], 'plot')
        ax1.text(-20, 10, r'$r_{\rm{event}} = 1$ Mpc')
        xmin = min(np.log10(min(tgw_amp_avg_SMBH[0])), np.log10(min(tgw_amp_avg_SMBH[1])), 
                   np.log10(min(tgw_amp_avg_IMBH[0])), np.log10(min(tgw_amp_avg_IMBH[1])))
        xmax = max(np.log10(min(tgw_amp_avg_SMBH[0])), np.log10(min(tgw_amp_avg_SMBH[1])),
                   np.log10(min(tgw_amp_avg_IMBH[0])), np.log10(min(tgw_amp_avg_IMBH[1])))
        ax1.set_xlim(1.02*xmin, 0.95*xmax)
        ax2.set_xlim(1.02*xmin, 0.95*xmax)
        ax1.legend()
        plt.savefig('figures/gravitational_waves/avg_GWstrain_histogram_ecc.pdf', dpi = 300, bbox_inches='tight')
        plt.clf()

        fig = plt.figure(figsize=(12.5, 6))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax_ = [ax1, ax2]
        for int_ in range(2):
            ax_[int_].scatter(tgw_frq_min_SMBH[int_], tgw_amp_min_SMBH[int_], color = 'black', edgecolors='black',  label = r'SMBH-IMBH')
            ax_[int_].scatter(tgw_frq_min_IMBH[int_], tgw_amp_min_IMBH[int_], color = 'purple', edgecolors='black',  label = r'IMBH-IMBH')
            ax_[int_].set_yscale('log')
            ax_[int_].set_xlabel(r'$f$ [Hz]')
            ax_[int_].set_ylabel(r'$h$')
            ax_[int_].set_title(str(integrator[int_]))
            plot_ini.tickers(ax_[int_], 'plot')
            ax_[int_].set_ylim(10**-25, 10**-16)
            ax_[int_].set_xscale('log')
            ax_[int_].set_yscale('log')
        ax1.legend()
        plt.savefig('figures/gravitational_waves/GW_freq_strain_maximise_diagram.pdf', dpi = 300, bbox_inches='tight')

        ecc_SMBH = 0
        ecc_IMBH = 0
        sem_SMBH = self.coll_radius(self.mass_SMBH[int_][k][0])
        sem_IMBH = self.coll_radius(self.mass_IMBH[int_][i][0])
        avg_amp_SMBH, avg_frq_SMBH = self.gw_amp_freq(sem_SMBH, ecc_SMBH, self.mass_SMBH[int_][i][0], self.mass_SMBH[int_][i][1])
        avg_amp_IMBH, avg_frq_IMBH = self.gw_amp_freq(sem_IMBH, ecc_IMBH, self.mass_IMBH[int_][i][0], self.mass_IMBH[int_][i][1])

        with open('figures/gravitational_waves/output/Binaries_redshift_freq_strain.txt', 'w') as file:
            for int_ in range(2):
                file.write('\nData for '+str(integrator[int_]))
                file.write('\nFor eccentric events @ z ~ 3.5, IMBH-IMBH avg. strain is: ' + str(np.mean(tgw_amp_avg_IMBH[int_]) * (1/7015)))
                file.write('\nFor eccentric events @ z ~ 3.5, IMBH-IMBH avg. freq is:   ' + str(np.mean(tgw_frq_avg_IMBH[int_]) * (1/7015)))
                file.write('\nFor eccentric events @ z ~ 3.5, IMBH-IMBH min. strain is: ' + str(np.mean(tgw_amp_min_IMBH[int_]) * (1/7015)))
                file.write('\nFor eccentric events @ z ~ 3.5, IMBH-IMBH min. freq is:   ' + str(np.mean(tgw_frq_min_IMBH[int_]) * (1/7015)))
                file.write('\nFive smallest GW frequencies for IMBH-IMBH:               ' + str(np.sort(tgw_frq_avg_IMBH[int_][:5]) + ' Myr'))
                file.write('\nFive smallest GW strains for IMBH-IMBH:                   ' + str(np.sort(tgw_amp_avg_IMBH[int_][:5]) + ' Myr'))
                file.write('\nFive smallest GW frequencies for SMBH-IMBH:               ' + str(np.sort(tgw_frq_avg_SMBH[int_][:5]) + ' Myr'))
                file.write('\nFive smallest GW strains for SMBH-IMBH:                   ' + str(np.sort(tgw_amp_avg_SMBH[int_][:5]) + ' Myr'))
                file.write('\nFor inspiral events @ z ~ 3.5,  IMBH-IMBH strain is:      ' + str(avg_amp_IMBH * (1/7015)))
                file.write('\nFor inspiral events @ z ~ 3.5,  IMBH-IMBH freq is:        ' + str(avg_frq_IMBH * (1/7015)))
                file.write('\nFor eccentric events @ z ~ 3.5, SMBH-IMBH avg. strain is: ' + str(np.mean(tgw_amp_avg_SMBH[int_]) * (1/7015)))
                file.write('\nFor eccentric events @ z ~ 3.5, IMBH-IMBH avg. freq is:   ' + str(np.mean(tgw_frq_avg_SMBH[int_]) * (1/7015)))
                file.write('\nFor eccentric events @ z ~ 3.5, SMBH-IMBH min. strain is: ' + str(np.mean(tgw_amp_min_SMBH[int_]) * (1/7015)))
                file.write('\nFor eccentric events @ z ~ 3.5, IMBH-IMBH min. freq is:   ' + str(np.mean(tgw_frq_min_SMBH[int_]) * (1/7015)))
                file.write('\nFor inspiral events @ z ~ 3.5,  SMBH-IMBH strain is:      ' + str(avg_amp_SMBH * (1/7015)))
                file.write('\nFor inspiral events @ z ~ 3.5,  SMBH-IMBH freq is:        ' + str(avg_frq_SMBH * (1/7015)))
                file.write('\n========================================================================')

    def IMBH_tgw_plotter(self):
        """
        Function which detects all binaries and plots their eccentricity vs. semi-major axis
        """

        plot_init = plotter_setup()

        tgw_IMBH_avg = [[ ], [ ]]
        tgw_IMBH_min = [[ ], [ ]]
        semi_IMBH_min = [[ ], [ ]]
        semi_IMBH_avg = [[ ], [ ]]
        semi_IMBH_conc = [[ ], [ ]]
        ecc_IMBH_min = [[ ], [ ]]
        ecc_IMBH_avg = [[ ], [ ]]
        ecc_IMBH_conc = [[ ], [ ]]
        
        for int_ in range(2):
            for i in range(len(self.semi_IMBH_avg[int_])):
                grav_IMBH_avg_timescale = self.gw_calc(self.semi_IMBH_avg[int_][i], self.ecc_IMBH_avg[int_][i], 
                                                       self.mass_IMBH[int_][i][0], self.mass_IMBH[int_][i][1]).value_in(units.Myr)
                grav_IMBH_min_timescale = self.gw_calc(self.semi_IMBH_min[int_][i], self.ecc_IMBH_min[int_][i], 
                                                       self.mass_IMBH[int_][i][0], self.mass_IMBH[int_][i][1]).value_in(units.Myr)
                tgw_IMBH_avg[int_].append(grav_IMBH_avg_timescale)
                tgw_IMBH_min[int_].append(grav_IMBH_min_timescale)

            semi_IMBH_avg_raw = np.asarray([i.value_in(units.parsec) for i in self.semi_IMBH_avg[int_]])
            semi_IMBH_min_raw = np.asarray([i.value_in(units.parsec) for i in self.semi_IMBH_min[int_]])
            ecc_IMBH_avg_raw = np.asarray(self.ecc_IMBH_avg[int_])
            ecc_IMBH_min_raw = np.asarray(self.ecc_IMBH_min[int_])
            tgw_IMBH_min_raw = np.asarray(tgw_IMBH_min[int_])

            ecc_IMBH_avg[int_] = np.asarray([np.log10(1-i) for i in ecc_IMBH_avg_raw[(ecc_IMBH_avg_raw < 1) & (semi_IMBH_avg_raw > 0)]])
            semi_IMBH_avg[int_] = np.asarray([np.log10(i) for i in semi_IMBH_avg_raw[(ecc_IMBH_avg_raw < 1) & (semi_IMBH_avg_raw > 0)]])
            ecc_IMBH_min[int_] = np.asarray([np.log10(1-i) for i in ecc_IMBH_min_raw[(ecc_IMBH_min_raw < 1) & (semi_IMBH_min_raw > 0)]])
            semi_IMBH_min[int_] = np.asarray([np.log10(i) for i in semi_IMBH_min_raw[(ecc_IMBH_min_raw < 1) & (semi_IMBH_min_raw > 0)]])
            tgw_IMBH_min[int_] = np.asarray(tgw_IMBH_min_raw[(ecc_IMBH_min_raw < 1) & (semi_IMBH_min_raw > 0)])
            ecc_IMBH_min[int_] = np.asarray(ecc_IMBH_min[int_][semi_IMBH_min[int_] > -5])
            tgw_IMBH_min[int_] = np.asarray(tgw_IMBH_min[int_][semi_IMBH_min[int_] > -5])
            semi_IMBH_min[int_] = np.asarray(semi_IMBH_min[int_][semi_IMBH_min[int_] > -5])

            ecc_IMBH_conc[int_] = np.concatenate((ecc_IMBH_avg[int_], ecc_IMBH_min[int_]), axis = None)
            semi_IMBH_conc[int_] = np.concatenate((semi_IMBH_avg[int_], semi_IMBH_min[int_]), axis = None)        
        
        ############### PLOTTING OF a vs. (1-e) FOR BIN ##############

        integrator = ['Hermite', 'GRX']
        
        xmin = min(np.nanmin(semi_IMBH_conc[0]), np.nanmin(semi_IMBH_conc[1]))
        xmax = max(np.nanmax(semi_IMBH_conc[0]), np.nanmax(semi_IMBH_conc[1]))
        ymin = min(np.nanmin(ecc_IMBH_conc[0]), np.nanmin(ecc_IMBH_conc[1]))
        cmin = min(np.nanmin(tgw_IMBH_min[0]), np.nanmin(tgw_IMBH_min[1]))
        cmax = max(np.nanmax(tgw_IMBH_min[0]), np.nanmax(tgw_IMBH_min[1]))
        norm_min = np.log10(cmin)
        norm_max = np.log10(cmax)
        normalise = plt.Normalize(norm_min, norm_max)

        for int_ in range(2):
            fig = plt.figure(figsize=(15, 6))
            ax1 = fig.add_subplot(121)
            ax2 = fig.add_subplot(122)
            bin2d_sim, xedg, xedg, image = ax2.hist2d(semi_IMBH_conc[int_], ecc_IMBH_conc[int_], bins=(75,75), 
                                                      range=([1.1*xmin, 1.1*xmax], [1.1*ymin, 0]))
            idx_merge = np.where(tgw_IMBH_min[int_] <= (self.tH))
            idx_stab = np.where(tgw_IMBH_min[int_] > (self.tH))
            bin2d_sim /= np.max(bin2d_sim)
            extent = [1.1*xmin, 1.1*xmax, 1.1*ymin, 0]
            contours = ax2.imshow(np.log10(bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
            for ax_ in [ax1, ax2]:
                ax_.set_ylabel(r'$\log_{10}(1-e)$')
                ax_.set_xlabel(r'$\log_{10} a$ [pc]')
                ax_.set_ylim(1.1*ymin, 0)
                ax_.set_xlim(-10, 0)
                
            plot_init.tickers(ax2, 'hist')
            plot_init.tickers(ax1, 'plot')
            colour_axes = ax1.scatter((semi_IMBH_min[int_][idx_merge]), ecc_IMBH_min[int_][idx_merge], 
                                       c = np.log10(tgw_IMBH_min[int_][idx_merge]), 
                                       norm = normalise, edgecolors='black',  marker = 'X')
            colour_axes = ax1.scatter((semi_IMBH_min[int_][idx_stab]), ecc_IMBH_min[int_][idx_stab], 
                                       c = np.log10(tgw_IMBH_min[int_][idx_stab]), 
                                       norm = normalise, edgecolors='black')
            self.forecast_interferometer(ax1, self.mass_IMBH[int_])
            plt.colorbar(colour_axes, ax = ax1, label = r'$\log_{10} \langle t_{GW}\rangle$ [Myr]')
            ax1.text(-9, -0.2, r'$f = 200$ Hz', verticalalignment = 'center', fontsize ='xx-small', rotation=self.text_angle-2, color = 'black')
            ax1.text(-5.3, -0.2, r'$f = 10^{-2}$ Hz', verticalalignment = 'center', fontsize ='xx-small', rotation=self.text_angle-2, color = 'black')
            plt.savefig('figures/gravitational_waves/ecc_semi_bins_IMBH_histogram'+str(integrator[int_])+'.pdf', dpi=300, bbox_inches='tight')


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
            idx_merge = np.where(tgw_SMBH_min[int_] <= (self.tH))
            idx_stab = np.where(tgw_SMBH_min[int_] > (self.tH))
            
            norm_min = np.log10(min(np.nanmin(tgw_SMBH_min[0]), np.nanmin(tgw_SMBH_min[0]),
                                    np.nanmin(tgw_SMBH_min[1]), np.nanmin(tgw_SMBH_min[1])))
            norm_max = np.log10(max(np.nanmax(tgw_SMBH_min[0]), np.nanmax(tgw_SMBH_min[0]),
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
            ax1.set_xlim(-10, 0)
            ax1.text(-9, -0.2, r'$f = 200$ Hz', verticalalignment = 'center', fontsize ='xx-small', rotation=self.text_angle-2, color = 'black')
            ax1.text(-5.3, -0.2, r'$f = 10^{-2}$ Hz', verticalalignment = 'center', fontsize ='xx-small', rotation=self.text_angle-2, color = 'black')
            plt.savefig('figures/gravitational_waves/ecc_semi_bins_SMBH_'+str(integrator[int_])+'_histogram.pdf', dpi=300, bbox_inches='tight')

        with open('figures/gravitational_waves/output/GWmerger_time.txt', 'w') as file:
            for int_ in range(2):
                file.write('\nData for '+str(integrator[int_]))
                file.write('\nAverage GW timescales for IMBH-SMBH:       ' + str(np.mean(tgw_SMBH_min[int_])) + ' Myr')
                file.write('\nFive smallest GW timescales for IMBH-SMBH: ' + str(np.sort(tgw_SMBH_min[int_])[:5]) + ' Myr')
                file.write('\n========================================================================')

        fig = plt.figure(figsize=(14, 6))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax1.set_title('Hermite')
        ax2.set_title('GRX')

        xmax = max(np.nanmax((self.close_enc[0])), np.nanmax((self.close_enc[1])))
        ymin = min(np.nanmin(np.log10(tgw_ratio[0])), np.nanmin(np.log10(tgw_ratio[0])))
        ymax = max(np.nanmax(np.log10(tgw_ratio[0])), np.nanmax(np.log10(tgw_ratio[0])))
        iter = -1
        for ax_ in [ax1, ax2]:
            iter += 1
            bin2d_sim, xedg, xedg, image = ax_.hist2d((self.close_enc[iter]), np.log10(tgw_ratio[iter]), 
                                                      bins=(50,50), range=([0, 1.1*xmax], [0.9*ymin, 1.1*ymax]))
            bin2d_sim /= np.max(bin2d_sim)
            extent = [0, 1.1*xmax, 0.9*ymin, 1.1*ymax]
            contours = ax_.imshow((bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
            ax_.set_xlabel(r'$N_{\rm{enc}}$')
            ax_.set_ylabel(r'$\log_{10}(t_{\rm{GW},f}/t_{\rm{GW},0})$')
            ax_.set_xlim(0, 1.1*xmax)
            plot_init.tickers(ax_, 'hist')
        plt.savefig('figures/gravitational_waves/Nenc_tgw_SMBH_evol_histogram.pdf', dpi=300, bbox_inches='tight')
        plt.clf()

    def transient_events(self):

        plot_init = plotter_setup()

        xdata_set1 = [self.time_flyby_nn[0], self.time_flyby_nn[1]]
        ydata_set1 = [self.freq_flyby_nn[0], self.freq_flyby_nn[1]]
        xdata_set2 = [self.time_flyby_t[0], self.time_flyby_t[1]]
        ydata_set2 = [self.freq_flyby_t[0], self.freq_flyby_t[1]]      

        fig = plt.figure(figsize=(16, 14))
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        integrator = ['Hermite', 'GRX']
        plot_title = ['Fly By', 'Tertiary']
        colours = ['red', 'blue']
        ax = [ax1, ax2]
        iter = -1
        for j in range(2):
            iter += 1
            plot_init.tickers(ax[j], 'plot')
            ax[j].set_yscale('log')
            ax[j].set_xlabel('Time [Myr]')
            ax[j].set_ylabel(r' $f_{\rm{GW}}$ [Hz]')
            for i in range(2):
                if j == 0:
                    ax[j].set_title(plot_title[j])
                    ax[j].scatter(xdata_set1[i], ydata_set1[i], edgecolors = 'black', 
                                  color = colours[i], label = integrator[i])
                else:
                    ax[j].set_title(plot_title[j])
                    ax[j].scatter(xdata_set2[i], ydata_set2[i], edgecolors = 'black', 
                                  color = colours[i])

        ax1.legend()
        plt.savefig('figures/gravitational_waves/events_time.pdf', dpi = 300, bbox_inches='tight')
