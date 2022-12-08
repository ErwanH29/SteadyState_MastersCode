from amuse.lab import *
from file_logistics import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import warnings
import pandas as pd
import statsmodels.api as sm

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
        warnings.filterwarnings("ignore", category=UserWarning) 
        self.H0 = 73.04 #Taken from arXiv:2112.04510 in km/s/Mpc
        self.tH = (self.H0*(3.2408*10**-20))**-1 * 1/(3600*365*24*10**6)
        self.integrator = ['Hermite', 'GRX']
        
    def new_data_extractor(self):
        """
        Script to extract data from recently simulated runs
        """

        print('!!!!!! WARNING THIS WILL TAKE A WHILE !!!!!!!')

        Hermite_data = glob.glob(os.path.join('/media/erwanh/Elements/Hermite/particle_trajectory/*'))
        GRX_data = glob.glob(os.path.join('/media/erwanh/Elements/GRX/particle_trajectory/*'))
        filename = [natsort.natsorted(Hermite_data), natsort.natsorted(GRX_data)] 
        for int_ in range(2):
            for file_ in range(len(filename[int_])):
                with open(filename[int_][file_], 'rb') as input_file:
                    print('Reading file', file_, ': ', filename[int_][file_])
                    data = pkl.load(input_file)
                    sim_time = np.shape(data)[1]-1

                    for parti_ in range(np.shape(data)[0]):
                        count = len(fnmatch.filter(os.listdir('data/tGW/'), '*.*'))
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
                        Nclose_indiv = 0
                        Nfbnn_indiv = 0
                        Nfbt_indiv = 0

                        semi_IMBH_avgv = -5
                        semi_IMBH_minv = -5
                        ecc_IMBH_avgv = -5
                        ecc_IMBH_minv = -5
                        time_IMBH_avgv = -5
                        time_IMBH_minv = -5
                        freq_NN_GW_indiv = [-5]
                        strain_NN_GW_indiv = [-5]
                        time_NN_GW_indiv = [-5]
                        SMBH_NN_event = [-5]
                        freq_t_GW_indiv = [-5]
                        strain_t_GW_indiv = [-5]
                        time_t_GW_indiv = [-5]
                        SMBH_t_event = [-5]

                        if parti_ != 0:
                            for col_ in range(np.shape(data)[1]-1):
                                sem_SMBH = data.iloc[parti_][col_][7][0]
                                if data.iloc[parti_][col_][8][0] < 1:
                                    semi_SMBH.append(sem_SMBH.value_in(units.pc))
                                    ecc_SMBH.append(data.iloc[parti_][col_][8][0])
                                    val = (data.iloc[parti_][col_][7][0].value_in(units.pc))**4 * (1-data.iloc[parti_][col_][8][0]**2)**3.5
                                    min_GW_SMBH.append(val)
                                
                                semi_major_nn = abs(data.iloc[parti_][col_][7][1])
                                semi_major_t = abs(data.iloc[parti_][col_][7][2])
                                ecc_nn = (data.iloc[parti_][col_][8][1])
                                ecc_t = (data.iloc[parti_][col_][8][2])

                                for part_ in range(np.shape(data)[0]):
                                    if data.iloc[part_][0][0] == data.iloc[parti_][col_][6][1]:
                                        mass2 = data.iloc[part_][0][1]
                                
                                strain_nn, freqGW_nn = self.gw_amp_freq(semi_major_nn, ecc_nn, mass1, mass2)
                                linex = data.iloc[parti_][col_][2][0] - data.iloc[0][col_][2][0]
                                liney = data.iloc[parti_][col_][2][1] - data.iloc[0][col_][2][1]
                                linez = data.iloc[parti_][col_][2][2] - data.iloc[0][col_][2][2]
                                dist_SMBH = (linex**2+liney**2+linez**2).sqrt()
                                dist_NN = data.iloc[parti_][col_][-1]

                                if freqGW_nn > 5*10**-10:
                                    freq_NN_GW_indiv.append(freqGW_nn)
                                    strain_NN_GW_indiv.append(strain_nn)
                                    time_NN_GW_indiv.append(10**-3 * col_)

                                    if dist_SMBH.value_in(units.pc) == dist_NN:
                                        Nfbnn_indiv += 1
                                        SMBH_NN_event.append(1)
                                    else:
                                        Nfbnn_indiv += 0.5
                                        SMBH_NN_event.append(-5)
                                        
                                tSMBH = False
                                if sem_SMBH == semi_major_t:
                                    mass2 = SMBH_sys_mass
                                    tSMBH = True

                                strain_t, freqGW_t = self.gw_amp_freq(semi_major_t, ecc_t, mass1, mass2)
                                if freqGW_t > 5*10**-10:
                                    sem_SMBH = data.iloc[parti_][col_][7][0]
                                    freq_t_GW_indiv.append(freqGW_t)
                                    strain_t_GW_indiv.append(strain_t)
                                    time_t_GW_indiv.append(10**-3 * col_)     

                                    if (tSMBH):
                                        Nfbt_indiv += 1
                                        SMBH_t_event.append(1)
                                    else:
                                        Nfbt_indiv += 0.5
                                        SMBH_t_event.append(-5)

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
                                    Nclose_indiv += 1

                            mass_IMBH = np.asarray(mass_IMBH)
                            IMBH_in = 0
                            mass_IMBH_binmin = 0 | units.MSun
                            if len(IMBH_key) > 1:
                                IMBH_sem = np.asarray(IMBH_sem)
                                IMBH_ecc = np.asarray(IMBH_ecc)
                                for rep_ in range(len(IMBH_key)):
                                    if IMBH_key[rep_] != IMBH_key[rep_-1] and rep_ != len(IMBH_key)-1:
                                        semi_IMBH_avgv = np.mean(IMBH_sem[IMBH_in:rep_+1])
                                        ecc_IMBH_avgv = np.mean(IMBH_ecc[IMBH_in:rep_+1])
                                        time_IMBH_avgv = np.mean(bin_time[IMBH_in:rep_+1])

                                        idx = np.nanargmin(min_GW_bin)
                                        semi_IMBH_minv = IMBH_sem[idx]
                                        ecc_IMBH_minv = IMBH_ecc[idx]
                                        time_IMBH_minv = bin_time[idx]
                                        mass_IMBH_binmin = mass_IMBH[idx][1] * (1 | units.MSun)

                                        IMBH_in = rep_

                                    if IMBH_key[rep_] != IMBH_key[rep_-1] and rep_ == len(IMBH_key)-1: #Otherwise blocks out at the last time step
                                        semi_IMBH_avgv = np.mean(IMBH_sem[IMBH_in:rep_+1])
                                        ecc_IMBH_avgv = np.mean(IMBH_ecc[IMBH_in:rep_+1])
                                        time_IMBH_avgv = np.mean(bin_time[IMBH_in:rep_+1])

                                        idx = np.nanargmin(min_GW_bin)
                                        semi_IMBH_minv = IMBH_sem[idx]
                                        ecc_IMBH_minv = IMBH_ecc[idx]
                                        time_IMBH_minv = bin_time[idx]
                                        mass_IMBH_binmin = mass_IMBH[idx][1] * (1 | units.MSun)

                                        IMBH_in = rep_

                            if len(IMBH_key) == 1:
                                semi_IMBH_avgv = IMBH_sem[0]
                                ecc_IMBH_avgv = IMBH_ecc[0]
                                time_IMBH_avgv = bin_time[0]

                                semi_IMBH_minv = IMBH_sem[0]
                                ecc_IMBH_minv = IMBH_ecc[0]
                                time_IMBH_minv = bin_time[0]
                                mass_IMBH_binmin = mass_IMBH[0][1] * (1 | units.MSun)

                            Ntot_indiv = Nfbnn_indiv + Nfbt_indiv

                            semi_SMBH = np.asarray(semi_SMBH)
                            ecc_SMBH = np.asarray(ecc_SMBH)
                            idx_SMBH = np.nanargmin(min_GW_SMBH)

                            path = '/home/erwanh/Desktop/SteadyStateBH/Data_Process/data/tGW/'
                            stab_tracker = pd.DataFrame()
                            df_stabtime = pd.Series({'Integrator': self.integrator[int_],
                                                     'Simulation Time': 10**3*sim_time,
                                                     'Population': 10*round(0.1*np.shape(data)[0]),
                                                     'Semi avg SMBH': np.mean(semi_SMBH),
                                                     'Semi min SMBH': semi_SMBH[idx_SMBH],
                                                     'Ecc avg SMBH': np.mean(ecc_SMBH),
                                                     'Ecc min SMBH': ecc_SMBH[idx_SMBH],
                                                     'Semi avg IMBH': semi_IMBH_avgv,
                                                     'Semi min IMBH': semi_IMBH_minv,
                                                     'Ecc avg IMBH': ecc_IMBH_avgv,
                                                     'Ecc min IMBH': ecc_IMBH_minv,
                                                     'Event avg time IMBH': time_IMBH_avgv,
                                                     'Event min time IMBH': time_IMBH_minv,
                                                     'mass SMBH': SMBH_sys_mass,
                                                     'mass IMBH1': mass1,
                                                     'mass IMBH2': mass_IMBH_binmin,
                                                     'No. Binary Events': Nfbnn_indiv,
                                                     'No. Tertiary Events': Nfbt_indiv,
                                                     'No. Total Events': Ntot_indiv,
                                                     'No. Close Encounter': Nclose_indiv,
                                                     'FlyBy Binary Frequencies': freq_NN_GW_indiv,
                                                     'Flyby Binary Strain': strain_NN_GW_indiv,
                                                     'Flyby Binary Time': time_NN_GW_indiv,
                                                     'Flyby SMBH Event': SMBH_NN_event,
                                                     'Flyby Tertiary Frequencies': freq_t_GW_indiv,
                                                     'Flyby Tertiary Strain': strain_t_GW_indiv,
                                                     'Flyby Tertiary Time': time_t_GW_indiv,
                                                     'Tertiary SMBH Event': SMBH_t_event})
                            stab_tracker = stab_tracker.append(df_stabtime, ignore_index = True)
                            stab_tracker.to_pickle(os.path.join(path, 'IMBH_'+str(self.integrator[int_])+'_tGW_data_indiv_parti_'+str(count)+'.pkl'))

    def combine_data(self):
        """
        Function which extracts ALL data and provides 
        the merger rate of the system.
        """

        self.sim_time =[[ ], [ ]]
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

        self.tot_sim_time = [[ ], [ ]]
        self.tot_events = [[ ], [ ]]
        self.fb_nn_events = [[ ], [ ]]
        self.fb_nn_SMBH = [[ ], [ ]]
        self.fb_t_events = [[ ], [ ]]
        self.fb_t_SMBH = [[ ], [ ]]

        self.event_rate =[[ ], [ ]]

        tGW_data = natsort.natsorted(glob.glob('data/tGW/*'))
        for file_ in range(len(tGW_data)):
            with open(tGW_data[file_], 'rb') as input_file:
                data_file = pkl.load(input_file)
                if data_file.iloc[0][0] == 'Hermite':
                    int_idx = 0
                else:
                    int_idx = 1

                self.sim_time[int_idx].append(data_file.iloc[0][1])
                self.pop[int_idx].append(int(data_file.iloc[0][2]))
                self.semi_SMBH_avg[int_idx].append(data_file.iloc[0][3])
                self.semi_SMBH_min[int_idx].append(data_file.iloc[0][4])
                self.ecc_SMBH_avg[int_idx].append(data_file.iloc[0][5])
                self.ecc_SMBH_min[int_idx].append(data_file.iloc[0][6])
                self.semi_IMBH_avg[int_idx].append(data_file.iloc[0][7])
                self.semi_IMBH_min[int_idx].append(data_file.iloc[0][8])
                self.ecc_IMBH_avg[int_idx].append(data_file.iloc[0][9])
                self.ecc_IMBH_min[int_idx].append(data_file.iloc[0][10])
                self.time_IMBH_avg[int_idx].append(data_file.iloc[0][11])
                self.time_IMBH_min[int_idx].append(data_file.iloc[0][12])
                self.mass_SMBH[int_idx].append([data_file.iloc[0][14], data_file.iloc[0][13]])
                self.mass_IMBH[int_idx].append([data_file.iloc[0][14], data_file.iloc[0][15]])
                self.fb_nn_events[int_idx].append(data_file.iloc[0][16])
                self.fb_t_events[int_idx].append(data_file.iloc[0][17])
                self.tot_events[int_idx].append(data_file.iloc[0][18])
                self.close_enc[int_idx].append(data_file.iloc[0][19])
                self.freq_flyby_nn[int_idx].append(data_file.iloc[0][20])
                self.strain_flyby_nn[int_idx].append(data_file.iloc[0][21])
                self.time_flyby_nn[int_idx].append(data_file.iloc[0][22])
                self.fb_nn_SMBH[int_idx].append(data_file.iloc[0][23])
                self.freq_flyby_t[int_idx].append(data_file.iloc[0][24])
                self.strain_flyby_t[int_idx].append(data_file.iloc[0][25])
                self.time_flyby_t[int_idx].append(data_file.iloc[0][26])
                self.fb_t_SMBH[int_idx].append(data_file.iloc[0][27])
                self.event_rate[int_idx].append(data_file.iloc[0][18]/data_file.iloc[0][1])

        for int_ in range(2):
            self.sim_time[int_] = np.asarray(self.sim_time[int_], dtype = 'object')
            self.pop[int_] = np.asarray(self.pop[int_], dtype = 'object')
            self.semi_SMBH_avg[int_] = np.asarray(self.semi_SMBH_avg[int_], dtype = 'object')
            self.semi_SMBH_min[int_] = np.asarray(self.semi_SMBH_min[int_], dtype = 'object')
            self.ecc_SMBH_avg[int_] = np.asarray(self.ecc_SMBH_avg[int_], dtype = 'object')
            self.ecc_SMBH_min[int_] = np.asarray(self.ecc_SMBH_min[int_], dtype = 'object')
            self.semi_IMBH_avg[int_] = np.asarray(self.semi_IMBH_avg[int_], dtype = 'object')
            self.semi_IMBH_min[int_] = np.asarray(self.semi_IMBH_min[int_], dtype = 'object')
            self.ecc_IMBH_avg[int_] = np.asarray(self.ecc_IMBH_avg[int_], dtype = 'object')
            self.ecc_IMBH_min[int_] = np.asarray(self.ecc_IMBH_min[int_], dtype = 'object')
            self.time_IMBH_avg[int_] = np.asarray(self.time_IMBH_avg[int_], dtype = 'object')
            self.time_IMBH_min[int_] = np.asarray(self.time_IMBH_min[int_], dtype = 'object')
            self.mass_SMBH[int_] = np.asarray(self.mass_SMBH[int_], dtype = 'object')
            self.mass_IMBH[int_] = np.asarray(self.mass_IMBH[int_], dtype = 'object')
            self.fb_nn_events[int_] = np.asarray(self.fb_nn_events[int_], dtype = 'object')
            self.fb_t_events[int_] = np.asarray(self.fb_t_events[int_], dtype = 'object')
            self.tot_events[int_] = np.asarray(self.tot_events[int_], dtype = 'object')
            self.close_enc[int_] = np.asarray(self.close_enc[int_], dtype = 'object')
            self.freq_flyby_nn[int_] = np.asarray(self.freq_flyby_nn[int_], dtype = 'object')
            self.strain_flyby_nn[int_] = np.asarray(self.strain_flyby_nn[int_], dtype = 'object')
            self.time_flyby_nn[int_] = np.asarray(self.time_flyby_nn[int_], dtype = 'object')
            self.freq_flyby_t[int_] = np.asarray(self.freq_flyby_t[int_], dtype = 'object')
            self.strain_flyby_t[int_] = np.asarray(self.strain_flyby_t[int_], dtype = 'object')
            self.time_flyby_t[int_] = np.asarray(self.time_flyby_t[int_], dtype = 'object')
            self.event_rate[int_] = np.asarray(self.event_rate[int_], dtype = 'object')
        
        with open('figures/gravitational_waves/output/event_rate.txt', 'w') as file:
                for int_ in range(2):
                    pop_idx = np.where(self.pop[int_] > 10)
                    #if len(pop_idx) > 0:
                    tot_events_t = self.tot_events[int_]#[pop_idx]
                    tot_events = tot_events_t[tot_events_t > 0]
                    fb_nn_events = self.fb_nn_events[int_]#[pop_idx]
                    fb_nn_events = fb_nn_events[tot_events_t > 0]
                    fb_t_events = self.fb_t_events[int_]#[pop_idx]
                    fb_t_events = fb_t_events[tot_events_t > 0]
                    sim_time = self.sim_time[int_]#[pop_idx]
                    sim_time = sim_time[tot_events_t > 0]

                    if len(tot_events) > 0:
                        file.write('\nData for '+str(self.integrator[int_]))
                        file.write('\nAverage total event rate per Myr            ' + str(np.mean(tot_events/sim_time * 10**6)))
                        file.write('\nAverage nearest neigh. event rate per Myr   ' + str(np.mean(fb_nn_events/sim_time * 10**6)))
                        file.write('\nAverage tertiary event rate per Myr         ' + str(np.mean(fb_t_events/sim_time * 10**6)))
                        file.write('\n========================================================================')

    def coll_radius(self, mass_arr):
        return 3 * (2*constants.G*mass_arr)/(constants.c**2)

    def forecast_interferometer(self, ax, mass_arr):
        """
        Function to plot the LISA and aLIGO frequency range in Ge a vs. (1-e) plots
        """
        
        ecc_range = np.linspace(0.0001, (1-10**-8), 50)

        self.LIGO_semimaj_max = self.GW_freq(ecc_range[1:], 10000 | units.Hz, mass_arr[0][0], mass_arr[0][1])
        self.LIGO_semimaj_min = self.GW_freq(ecc_range[1:], 10 | units.Hz, mass_arr[0][0], mass_arr[0][1])
        self.LIGO_semimaj = self.GW_freq(ecc_range[1:], 200 | units.Hz, mass_arr[0][0], mass_arr[0][1])

        self.LISA_semimaj_max = self.GW_freq(ecc_range[1:], 1 | units.Hz, mass_arr[0][0], mass_arr[0][1])
        self.LISA_semimaj_min = self.GW_freq(ecc_range[1:], 1e-4 | units.Hz, mass_arr[0][0], mass_arr[0][1])
        self.LISA_semimaj = self.GW_freq(ecc_range[1:], 1e-2 | units.Hz, mass_arr[0][0], mass_arr[0][1])

        ecc_range = [np.log(1-i) for i in ecc_range[1:]]
        self.text_angle = np.degrees(np.arctan((ecc_range[30]-ecc_range[20])/(self.LIGO_semimaj[30]-self.LIGO_semimaj[20]))) - 20

        ax.plot(self.LIGO_semimaj_min, ecc_range, linestyle = ':', color = 'black')
        ax.plot(self.LIGO_semimaj, ecc_range, linestyle = '-.', color = 'black')
        ax.plot(self.LIGO_semimaj_max, ecc_range, linestyle = ':', color = 'black')
        ax.fill_between(np.append(self.LIGO_semimaj_min, self.LIGO_semimaj_max[::-1]), 
                        np.append(ecc_range[:], ecc_range[::-1]), alpha = 0.2, color = 'blue')
        ax.plot(self.LISA_semimaj_min, ecc_range, linestyle = ':', color = 'black')
        ax.plot(self.LISA_semimaj, ecc_range, linestyle = '-.', color = 'black')
        ax.plot(self.LISA_semimaj_max, ecc_range, linestyle = ':', color = 'black')
        ax.fill_between(np.append(self.LISA_semimaj_min, self.LISA_semimaj_max[::-1]), 
                        np.append(ecc_range[:], ecc_range[::-1]), alpha = 0.2, color = 'red')

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

    def LISA_sensitivity(self, freq):
        """
        Function to plot the LISA sensitivity curve.
        Data are taken from arXiv:1803.01944
        
        Inputs:
        freq:    The frequency range which LISA occupies"""
        length = 2.5*10**9
        fstar  = 19.09*10**-3
        
        alpha = 0.138
        beta = -221
        kappa = 521
        gamma = 1680
        fk = 0.00113
        A = 9 * 10**-45

        R = 3./20./(1. + 6./10.*(freq/fstar)**2)*2

        p_oms = (1.5*10**-11)**2 * (1+((2*10**-3)/freq)**4)
        p_acc = (3*10**-15)**2 * (1+((0.4*10**-3)/freq)**2)*(1+(freq/(8*10**-3))**4)
        sens_const = A * freq**(-(7.3))*np.exp(-freq**alpha+beta*freq*np.sin(kappa*freq))*(1+np.tanh(gamma*(fk-freq)))
        Pn = (p_oms + 2.*(1. + np.cos(freq/fstar)**2)*p_acc/(2.*np.pi*freq)**4)/length**2

        return Pn/R + sens_const

    def scatter_hist(self, x, y, ax, ax_histf, ax_histh, len_arr):
        """
        Function to plot the frequency/strain histogram along its scatter plot
        
        Inputs:
        x:        The strain array
        y:        The frequency array
        ax:       axis where the scatter plot is located
        ax_histf: axis where the strain histogram is placed
        ax_histh: axis where the frequency histogram is placed
        len_arr:  Length of SMBH data array
        """
        x_temp = np.linspace(10**-5, 1, 10**3)
        PTA_freq = np.linspace(-9, -6)
        # the scatter plot:
        ax.plot(np.log10(x_temp), np.log10(np.sqrt(x_temp*self.LISA_sensitivity(x_temp))), color = 'black')
        ax.fill_between(PTA_freq, -26, -13, alpha = 0.35, color = 'red')
        ax.text(-2.3, -20.5, 'Lisa Sensitivity', fontsize ='small', rotation = 276)
        ax.text(-6, -24, 'PTA Window', fontsize ='small', rotation = 270)
        ax.scatter(np.log10(x[:len_arr]), np.log10(y[:len_arr]), label = 'SMBH-IMBH', color = 'black')
        ax.scatter(np.log10(x[len_arr:]), np.log10(y[len_arr:]), label = 'IMBH-IMBH', color = 'purple', edgecolors = 'black')
        
        ax_histf.tick_params(axis="x", labelbottom=False)
        ax_histh.tick_params(axis="y", labelleft=False)

        kdef_SMBH = sm.nonparametric.KDEUnivariate(np.log10(x[:len_arr]))
        kdef_SMBH.fit()
        kdef_SMBH.density /= max(kdef_SMBH.density)
        kdef_IMBH = sm.nonparametric.KDEUnivariate(np.log10(x[len_arr:]))
        kdef_IMBH.fit()
        kdef_IMBH.density /= max(kdef_IMBH.density)

        ax_histf.plot(kdef_SMBH.support, kdef_SMBH.density, color = 'black')
        ax_histf.fill_between(kdef_SMBH.support, kdef_SMBH.density, alpha = 0.35, color = 'black')
        ax_histf.plot(kdef_IMBH.support, kdef_IMBH.density, color = 'purple')
        ax_histf.fill_between(kdef_IMBH.support, kdef_IMBH.density, alpha = 0.35, color = 'purple')
        ax_histf.set_ylim(0, 1.05)

        kdeh_SMBH = sm.nonparametric.KDEUnivariate(np.log10(y[:len_arr]))
        kdeh_SMBH.fit()
        kdeh_SMBH.density /= max(kdeh_SMBH.density)
        kdeh_IMBH = sm.nonparametric.KDEUnivariate(np.log10(y[len_arr:]))
        kdeh_IMBH.fit()
        kdeh_IMBH.density /= max(kdeh_IMBH.density)

        ax_histh.plot(kdeh_SMBH.density, kdeh_SMBH.support, color = 'black')
        ax_histh.fill_between(kdeh_SMBH.density, kdeh_SMBH.support, alpha = 0.35, color = 'black')
        ax_histh.plot(kdeh_IMBH.density, kdeh_IMBH.support, color = 'purple')
        ax_histh.fill_between(kdeh_IMBH.density, kdeh_IMBH.support, alpha = 0.35, color = 'purple')
        ax_histh.set_xlim(0, 1.05)
        
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
        
        for int_ in range(1):
            semi_IMBH_avg_temp = []
            semi_IMBH_min_temp = []
            ecc_IMBH_avg_temp = []
            ecc_IMBH_min_temp = []
            for i in range(len(self.semi_IMBH_avg[int_])):
                if self.ecc_IMBH_avg[int_][i] > 0:
                    semi_avg = self.semi_IMBH_avg[int_][i] 
                    semi_min = self.semi_IMBH_min[int_][i]
                    semi_IMBH_avg_temp.append(semi_avg)
                    semi_IMBH_min_temp.append(semi_min)

                    ecc_avg = self.ecc_IMBH_avg[int_][i]
                    ecc_min = self.ecc_IMBH_min[int_][i]
                    ecc_IMBH_avg_temp.append(ecc_avg)
                    ecc_IMBH_min_temp.append(ecc_min)

                    grav_IMBH_avg_timescale = self.gw_calc(semi_avg * (1 | units.pc), ecc_avg, self.mass_IMBH[int_][i][0], 
                                                           self.mass_IMBH[int_][i][1]).value_in(units.Myr)
                    grav_IMBH_min_timescale = self.gw_calc(semi_min * (1 | units.pc), ecc_min, self.mass_IMBH[int_][i][0], 
                                                           self.mass_IMBH[int_][i][1]).value_in(units.Myr)
                    tgw_IMBH_avg[int_].append(grav_IMBH_avg_timescale)
                    tgw_IMBH_min[int_].append(grav_IMBH_min_timescale)

            semi_IMBH_avg_raw = np.asarray([i for i in semi_IMBH_avg_temp])
            semi_IMBH_min_raw = np.asarray([i for i in semi_IMBH_min_temp])
            ecc_IMBH_avg_raw = np.asarray(ecc_IMBH_avg_temp)
            ecc_IMBH_min_raw = np.asarray(ecc_IMBH_min_temp)
            tgw_IMBH_min_raw = np.asarray(tgw_IMBH_min[int_])
            tgw_IMBH_avg_raw = np.asarray(tgw_IMBH_avg[int_])

            ecc_IMBH_avg[int_] = np.asarray([np.log10(1-i) for i in ecc_IMBH_avg_raw[(ecc_IMBH_avg_raw < 1) & (semi_IMBH_avg_raw > 0)]])
            semi_IMBH_avg[int_] = np.asarray([np.log10(i) for i in semi_IMBH_avg_raw[(ecc_IMBH_avg_raw < 1) & (semi_IMBH_avg_raw > 0)]])
            ecc_IMBH_min[int_] = np.asarray([np.log10(1-i) for i in ecc_IMBH_min_raw[(ecc_IMBH_min_raw < 1) & (semi_IMBH_min_raw > 0)]])
            semi_IMBH_min[int_] = np.asarray([np.log10(i) for i in semi_IMBH_min_raw[(ecc_IMBH_min_raw < 1) & (semi_IMBH_min_raw > 0)]])
            tgw_IMBH_min[int_] = np.asarray(tgw_IMBH_min_raw[(ecc_IMBH_min_raw < 1) & (semi_IMBH_min_raw > 0)])
            tgw_IMBH_avg[int_] = np.asarray(tgw_IMBH_avg_raw[(ecc_IMBH_min_raw < 1) & (semi_IMBH_min_raw > 0)])
            ecc_IMBH_min[int_] = np.asarray(ecc_IMBH_min[int_][semi_IMBH_min[int_] > -5])
            tgw_IMBH_min[int_] = np.asarray(tgw_IMBH_min[int_][semi_IMBH_min[int_] > -5])
            tgw_IMBH_avg[int_] = np.asarray(tgw_IMBH_avg[int_][semi_IMBH_min[int_] > -5])
            semi_IMBH_min[int_] = np.asarray(semi_IMBH_min[int_][semi_IMBH_min[int_] > -5])

            ecc_IMBH_conc[int_] = np.concatenate((ecc_IMBH_avg[int_], ecc_IMBH_min[int_]), axis = None)
            semi_IMBH_conc[int_] = np.concatenate((semi_IMBH_avg[int_], semi_IMBH_min[int_]), axis = None)        
        
        ############### PLOTTING OF a vs. (1-e) FOR BIN ##############

        """xmin = min(np.nanmin(semi_IMBH_conc[0]), np.nanmin(semi_IMBH_conc[1]))
        xmax = max(np.nanmax(semi_IMBH_conc[0]), np.nanmax(semi_IMBH_conc[1]))
        ymin = min(np.nanmin(ecc_IMBH_conc[0]), np.nanmin(ecc_IMBH_conc[1]))
        cmin = min(np.nanmin(tgw_IMBH_min[0]), np.nanmin(tgw_IMBH_min[1]))
        cmax = max(np.nanmax(tgw_IMBH_min[0]), np.nanmax(tgw_IMBH_min[1]))"""

        xmin = min((semi_IMBH_conc[0]))
        xmax = max((semi_IMBH_conc[0]))
        ymin = min((ecc_IMBH_conc[0]))
        cmin = min((tgw_IMBH_min[0]))
        cmax = max((tgw_IMBH_min[0]))
        norm_min = np.log10(cmin)
        norm_max = np.log10(cmax)
        normalise = plt.Normalize(norm_min, norm_max)

        for int_ in range(1):
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
                ax_.set_xlim(-11, 1)
                
            plot_init.tickers(ax2, 'hist')
            plot_init.tickers(ax1, 'plot')
            colour_axes = ax1.scatter((semi_IMBH_min[int_][idx_merge]), ecc_IMBH_min[int_][idx_merge], 
                                       c = np.log10(tgw_IMBH_min[int_][idx_merge]), 
                                       norm = normalise, edgecolors='black',  marker = 'X')
            colour_axes = ax1.scatter((semi_IMBH_min[int_][idx_stab]), ecc_IMBH_min[int_][idx_stab], 
                                       c = np.log10(tgw_IMBH_min[int_][idx_stab]), 
                                       norm = normalise, edgecolors='black')
            self.forecast_interferometer(ax1, self.mass_IMBH[int_])
            plt.colorbar(colour_axes, ax = ax1, label = r'$t_{\rm{min, GW}}$ [Myr]')
            ax1.text(-10.5, -0.4, r'aLIGO ($f_{\rm{peak}} = 200$ Hz)', verticalalignment = 'center', fontsize ='small', rotation=self.text_angle, color = 'black')
            ax1.text(-7.65, -0.4, r'LISA ($f_{\rm{peak}} = 10^{-2}$ Hz)', verticalalignment = 'center', fontsize ='small', rotation=self.text_angle, color = 'black')
            plt.savefig('figures/gravitational_waves/ecc_semi_bins_IMBH_histogram'+str(self.integrator[int_])+'.pdf', dpi=300, bbox_inches='tight')
            plt.clf()

        with open('figures/gravitational_waves/output/GW_IMBHmerger_time.txt', 'w') as file:
            for int_ in range(2):
                file.write('\nData for '+str(self.integrator[int_]))
                file.write('\nAverage avg. GW timescales for IMBH-IMBH:           ' + str(np.mean(tgw_IMBH_avg[int_])) + ' Myr')
                file.write('\nFive smallest avg. GW timescales for IMBH-IMBH: ' + str(np.sort(tgw_IMBH_avg[int_])[:5]) + ' Myr')
                file.write('\nFive smallest GW timescales for IMBH-IMBH:     ' + str(np.sort(tgw_IMBH_min[int_])[:5]) + ' Myr')
                file.write('\n========================================================================')

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
        tgw_SMBH_min = [[ ], [ ]]
        tgw_SMBH_avg = [[ ], [ ]]
        ecc_SMBH_conc = [[ ], [ ]]
        sem_SMBH_conc = [[ ], [ ]]
        ecc_SMBH_min = [[ ], [ ]]
        ecc_SMBH_avg = [[ ], [ ]]
        semi_SMBH_min = [[ ], [ ]]
        semi_SMBH_avg = [[ ], [ ]]

        ymin_ecc_sem = [ ]
        xmin_ecc_sem = [ ] ; xmax_ecc_sem = [ ]

        for int_ in range(1):
            for i in range(len(self.semi_SMBH_avg[int_])):
                grav_avg_time = self.gw_calc(self.semi_SMBH_avg[int_][i] * (1 | units.pc), self.ecc_SMBH_avg[int_][i], 
                                             self.mass_SMBH[int_][i][0], self.mass_SMBH[int_][i][1]).value_in(units.Myr)
                tgw_SMBH_avg[int_].append(grav_avg_time)
                grav_min_time = self.gw_calc(self.semi_SMBH_min[int_][i] * (1 | units.pc), self.ecc_SMBH_min[int_][i], 
                                             self.mass_SMBH[int_][i][0], self.mass_SMBH[int_][i][1]).value_in(units.Myr)
                tgw_SMBH_min[int_].append(grav_min_time)

            tgw_SMBH_min_raw_arr = np.asarray(tgw_SMBH_min[int_])
            tgw_SMBH_avg_raw_arr = np.asarray(tgw_SMBH_avg[int_])
            ecc_SMBH_avg_raw_arr = np.asarray(self.ecc_SMBH_avg[int_])
            ecc_SMBH_min_raw_arr = np.asarray(self.ecc_SMBH_min[int_])
            semi_SMBH_avg_raw_arr = np.asarray([i for i in self.semi_SMBH_avg[int_]])
            semi_SMBH_min_raw_arr = np.asarray([i for i in self.semi_SMBH_min[int_]])

            ecc_SMBH_avg[int_] = np.asarray([np.log10(1-i) for i in ecc_SMBH_avg_raw_arr[(ecc_SMBH_avg_raw_arr< 1) & (semi_SMBH_avg_raw_arr > 0)]])
            semi_SMBH_avg[int_] = np.asarray([np.log(i) for i in semi_SMBH_avg_raw_arr[(ecc_SMBH_avg_raw_arr < 1) & (semi_SMBH_avg_raw_arr > 0)]])
            ecc_SMBH_min[int_] = np.asarray([np.log10(1-i) for i in ecc_SMBH_min_raw_arr[(ecc_SMBH_min_raw_arr < 1) & (semi_SMBH_min_raw_arr > 0)]])
            semi_SMBH_min[int_] = np.asarray([np.log10(i) for i in semi_SMBH_min_raw_arr[(ecc_SMBH_min_raw_arr < 1) & (semi_SMBH_min_raw_arr > 0)]])
            tgw_SMBH_min[int_] = np.asarray(tgw_SMBH_min_raw_arr[(ecc_SMBH_min_raw_arr < 1) & (semi_SMBH_min_raw_arr > 0)])
            tgw_SMBH_avg[int_] = np.asarray(tgw_SMBH_avg_raw_arr[(ecc_SMBH_min_raw_arr < 1) & (semi_SMBH_min_raw_arr > 0)])

            ecc_SMBH_conc[int_] = np.concatenate((ecc_SMBH_avg[int_], ecc_SMBH_min[int_]), axis = None)
            sem_SMBH_conc[int_] = np.concatenate((semi_SMBH_avg[int_], semi_SMBH_min[int_]), axis = None)

            xmin_ecc_sem.append(np.nanmin(semi_SMBH_min[int_]))
            xmax_ecc_sem.append(np.nanmax(semi_SMBH_min[int_]))
            ymin_ecc_sem.append(np.nanmin(ecc_SMBH_min[int_]))

        ############## PLOTTING OF a vs. (1-e) FOR SMBH ##############

        fig, ax = plt.subplots(figsize=(12,6), nrows=1, ncols=2)
        gs = gridspec.GridSpec(1, 2)
        ax1 = plt.subplot(gs[0,0])
        ax2 = plt.subplot(gs[0,1])
        ax1.set_title('Hermite')
        ax2.set_title('GRX')        
        for int_ in range(1):
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
            idx_merge = np.where(tgw_SMBH_avg[int_] <= (self.tH))
            idx_stab = np.where(tgw_SMBH_avg[int_] > (self.tH))
            """
            norm_min = np.log10(min(np.nanmin(tgw_SMBH_avg[0]), np.nanmin(tgw_SMBH_avg[1])))
            norm_max = np.log10(max(np.nanmax(tgw_SMBH_avg[0]), np.nanmax(tgw_SMBH_avg[1])))"""
            norm_min = np.log10((np.nanmin(tgw_SMBH_avg[0])))
            norm_max = np.log10(np.nanmax(tgw_SMBH_avg[0]))
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
            ax1.set_xlim(-10, 1.1*max(xmax_ecc_sem))
            ax1.text(-8.95, -2, r'aLIGO ($f_{\rm{peak}} = 200$ Hz)', verticalalignment = 'center', fontsize ='small', rotation=self.text_angle+9.05, color = 'black')
            ax1.text(-6.05, -2, r'LISA ($f_{\rm{peak}} = 10^{-2}$ Hz)', verticalalignment = 'center', fontsize ='small', rotation=self.text_angle+9.05, color = 'black')
            plt.savefig('figures/gravitational_waves/ecc_semi_bins_SMBH_'+str(self.integrator[int_])+'_histogram.pdf', dpi=300, bbox_inches='tight')

        with open('figures/gravitational_waves/output/GW_SMBHmerger_time.txt', 'w') as file:
            for int_ in range(2):
                file.write('\nData for '+str(self.integrator[int_]))
                file.write('\nAverage avg. GW timescales for IMBH-SMBH:           ' + str(np.mean(tgw_SMBH_avg[int_])) + ' Myr')
                file.write('\nFive smallest avg. GW timescales for IMBH-SMBH: ' + str(np.sort(tgw_SMBH_avg[int_])[:5]) + ' Myr')
                file.write('\nFive smallest GW timescales for IMBH-SMBH:     ' + str(np.sort(tgw_SMBH_min[int_])[:5]) + ' Myr')
                file.write('\n========================================================================')

        fig = plt.figure(figsize=(14, 6))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax1.set_title('Hermite')
        ax2.set_title('GRX')
        """
        xmax = max(np.nanmax((self.close_enc[0])), np.nanmax((self.close_enc[1])))

        ymin = min(np.nanmin(np.log10(tgw_SMBH_min[0])), np.nanmin(np.log10(tgw_SMBH_min[1])))
        ymax = max(np.nanmax(np.log10(tgw_SMBH_min[0])), np.nanmax(np.log10(tgw_SMBH_min[1])))"""


        xmax = np.log10(max(((self.close_enc[0]))))
        ymin = (np.nanmin(np.log10(tgw_SMBH_min[0])))
        ymax = (np.nanmax(np.log10(tgw_SMBH_min[0])))
        iter = -1
        for ax_ in [ax1, ax2]:
            iter += 1
            x_vals = [np.log(i) for i in self.close_enc[iter]]
            bin2d_sim, xedg, xedg, image = ax_.hist2d(x_vals, np.log10(tgw_SMBH_min[iter]), 
                                                      bins=(40,40), range=([0, 1.1*xmax], [0.9*ymin, 1.1*ymax]))
            bin2d_sim /= np.max(bin2d_sim)
            extent = [0, 1.1*xmax, 0.9*ymin, 1.1*ymax]
            contours = ax_.imshow((bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
            ax_.set_xlabel(r'$N_{\rm{enc}}$')
            ax_.set_ylabel(r'$\log_{10}t_{\rm{GW, min}}$')
            ax_.set_xlim(0, 1.1*xmax)
            plot_init.tickers(ax_, 'hist')
        plt.savefig('figures/gravitational_waves/Nenc_tgw_SMBH_evol_histogram.pdf', dpi=300, bbox_inches='tight')
        plt.clf()

    def strain_freq_plotter(self):
        """
        Function which plots the amplitude histogram for all sets of particle 
        using Gultekin et al. 2005 eqn 1.
        """

        plot_ini = plotter_setup()

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
                if self.ecc_IMBH_avg[int_][i] > 0:
                    ecc_IMBH = self.ecc_IMBH_avg[int_][i]
                    sem_IMBH = self.semi_IMBH_avg[int_][i] * (1 | units.parsec)
                    avg_amp, avg_freq = self.gw_amp_freq(sem_IMBH, ecc_IMBH, self.mass_IMBH[int_][i][0], self.mass_IMBH[int_][i][1])
                    tgw_amp_avg_IMBH[int_].append(avg_amp)
                    tgw_frq_avg_IMBH[int_].append(avg_freq)

                    ecc_IMBHm = self.ecc_IMBH_min[int_][i]
                    sem_IMBHm = self.semi_IMBH_min[int_][i] * (1 | units.parsec)
                    min_amp, min_freq = self.gw_amp_freq(sem_IMBHm, ecc_IMBHm, self.mass_IMBH[int_][i][0], self.mass_IMBH[int_][i][1])
                    tgw_amp_min_IMBH[int_].append(min_amp)
                    tgw_frq_min_IMBH[int_].append(min_freq)

            for k in range(len(self.semi_SMBH_min[int_])):
                if self.semi_SMBH_min[int_][k] > 0:
                    ecc_SMBH = self.ecc_SMBH_avg[int_][k]
                    sem_SMBH = self.semi_SMBH_avg[int_][k] * (1 | units.parsec)
                    avg_amp, avg_freq = self.gw_amp_freq(sem_SMBH, ecc_SMBH, self.mass_SMBH[int_][k][0], self.mass_SMBH[int_][k][1])
                    tgw_amp_avg_SMBH[int_].append(avg_amp)
                    tgw_frq_avg_SMBH[int_].append(avg_freq)

                    ecc_SMBHm = self.ecc_SMBH_min[int_][k]
                    sem_SMBHm = self.semi_SMBH_min[int_][k] * (1 | units.parsec)
                    min_amp, min_freq = self.gw_amp_freq(sem_SMBHm, ecc_SMBHm, self.mass_SMBH[int_][k][0], self.mass_SMBH[int_][k][1])
                    tgw_amp_min_SMBH[int_].append(min_amp)
                    tgw_frq_min_SMBH[int_].append(min_freq)

            tgw_amp_avg_IMBH[int_] = np.asarray(tgw_amp_avg_IMBH[int_])
            tgw_frq_avg_IMBH[int_] = np.asarray(tgw_frq_avg_IMBH[int_])
            tgw_amp_avg_SMBH[int_] = np.asarray(tgw_amp_avg_SMBH[int_])
            tgw_frq_avg_SMBH[int_] = np.asarray(tgw_frq_avg_SMBH[int_])

        fig = plt.figure(figsize=(8, 6))
        gs = fig.add_gridspec(2, 2,  width_ratios=(4, 2), height_ratios=(2, 4),
                              left=0.1, right=0.9, bottom=0.1, top=0.9,
                              wspace=0.05, hspace=0.05)
        ax = fig.add_subplot(gs[1, 0])
        ax1 = fig.add_subplot(gs[0, 0], sharex=ax)
        ax2 = fig.add_subplot(gs[1, 1], sharey=ax)
        for int_ in range(1):
            length_SMBH = len(tgw_frq_min_SMBH[int_])
            self.scatter_hist(np.concatenate((tgw_frq_min_SMBH[int_], tgw_frq_min_IMBH[int_])), 
                              np.concatenate((tgw_amp_min_SMBH[int_], tgw_amp_min_IMBH[int_])),
                              ax, ax1, ax2, length_SMBH)
            ax.set_xlabel(r'$\log_{10}f$ [Hz]')
            ax.set_ylabel(r'$\log_{10}h$')
            ax1.set_title(str(self.integrator[int_]))
            plot_ini.tickers(ax1, 'plot')
            plot_ini.tickers(ax2, 'plot')
            plot_ini.tickers(ax, 'plot')
            ax.set_ylim(-26, -13.2)
            ax.set_xlim(-12.5, 1.5)
            ax.legend()
            plt.savefig('figures/gravitational_waves/'+str(self.integrator[int_])+'GW_freq_strain_maximise_diagram.pdf', dpi = 300, bbox_inches='tight')
            plt.clf()

        ecc_SMBH = 0
        ecc_IMBH = 0
        SMBH_mass = self.mass_SMBH[0][0][0]
        IMBH_mass = self.mass_SMBH[0][1][0]
        sem_SMBH = self.coll_radius(SMBH_mass)
        sem_IMBH = self.coll_radius(IMBH_mass)
        avg_amp_SMBH, avg_frq_SMBH = self.gw_amp_freq(sem_SMBH, ecc_SMBH, SMBH_mass, SMBH_mass)
        avg_amp_IMBH, avg_frq_IMBH = self.gw_amp_freq(sem_IMBH, ecc_IMBH, IMBH_mass, IMBH_mass)

        with open('figures/gravitational_waves/output/Binaries_redshift_freq_strain.txt', 'w') as file:
            for int_ in range(2):
                file.write('\nData for '+str(self.integrator[int_]))
                file.write('\nFor eccentric events @ z ~ 3.5, SMBH-IMBH avg. strain is: ' + str(np.mean(tgw_amp_avg_SMBH[int_]) * (1/7015)))
                file.write('\nFor eccentric events @ z ~ 3.5, IMBH-IMBH avg. strain is: ' + str(np.mean(tgw_amp_avg_IMBH[int_]) * (1/7015)))
                file.write('\nFor eccentric events @ z ~ 3.5, SMBH-IMBH avg. freq is:   ' + str(np.mean(tgw_frq_avg_SMBH[int_]) * (1/7015)))
                file.write('\nFor eccentric events @ z ~ 3.5, IMBH-IMBH avg. freq is:   ' + str(np.mean(tgw_frq_avg_IMBH[int_]) * (1/7015)))
                file.write('\nFor eccentric events @ z ~ 3.5, SMBH-IMBH min. strain is: ' + str(np.mean(tgw_amp_min_SMBH[int_]) * (1/7015)))
                file.write('\nFor eccentric events @ z ~ 3.5, IMBH-IMBH min. strain is: ' + str(np.mean(tgw_amp_min_IMBH[int_]) * (1/7015)))
                file.write('\nFor eccentric events @ z ~ 3.5, SMBH-IMBH min. freq is:   ' + str(np.mean(tgw_frq_min_SMBH[int_]) * (1/7015)))
                file.write('\nFor eccentric events @ z ~ 3.5, IMBH-IMBH min. freq is:   ' + str(np.mean(tgw_frq_min_IMBH[int_]) * (1/7015)))
                file.write('\nThree largest GW freq. for IMBH-IMBH:                     ' + str(np.sort(tgw_frq_avg_IMBH[int_])[-3:]) + ' Myr')
                file.write('\nThree largest GW strains for IMBH-IMBH:                   ' + str(np.sort(tgw_amp_avg_IMBH[int_])[-3:]) + ' Myr')
                file.write('\nThree largest GW freq. for SMBH-IMBH:                     ' + str(np.sort(tgw_frq_avg_SMBH[int_])[-3:]) + ' Myr')
                file.write('\nThree largest GW strains for SMBH-IMBH:                   ' + str(np.sort(tgw_amp_avg_SMBH[int_])[-3:]) + ' Myr')
                file.write('\nFor inspiral events @ z ~ 3.5,  IMBH-IMBH strain is:      ' + str(avg_amp_IMBH * (1/7015)))
                file.write('\nFor inspiral events @ z ~ 3.5,  IMBH-IMBH freq is:        ' + str(avg_frq_IMBH * (1/7015)))
                file.write('\nFor inspiral events @ z ~ 3.5,  SMBH-IMBH strain is:      ' + str(avg_amp_SMBH * (1/7015)))
                file.write('\nFor inspiral events @ z ~ 3.5,  SMBH-IMBH freq is:        ' + str(avg_frq_SMBH * (1/7015)))
                file.write('\n========================================================================')

    def transient_events(self):

        plot_init = plotter_setup() 

        xdata_set1 = [[ ], [ ]]
        ydata_set1 = [[ ], [ ]]
        xdata_set2 = [[ ], [ ]]
        ydata_set2 = [[ ], [ ]]

        fig = plt.figure(figsize=(16, 14))
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        plot_title = ['Fly By', 'Tertiary']
        colours = ['red', 'blue']
        ax = [ax1, ax2]
        iter = -1
        for int_ in range(1):
            time_nn_arr = [ ]
            freq_nn_arr = [ ]
            marker_nn_arr = [ ]
            time_t_arr = [ ]
            freq_t_arr = [ ]
            marker_t_arr = [ ]

            iter += 1
            for i in range(len(self.time_flyby_nn[int_])):
                if self.pop[int_][i] == 30:
                    for k in range(len(self.time_flyby_nn[int_][i])):
                        if self.freq_flyby_nn[int_][i][k] > 0:
                            time_nn_arr.append(self.time_flyby_nn[int_][i][k])
                            freq_nn_arr.append(self.freq_flyby_nn[int_][i][k])
                            marker_nn_arr.append(self.fb_nn_SMBH[int_][i][k])

            for i in range(len(self.time_flyby_t[int_])):
                if self.pop[int_][i] == 30:
                    for k in range(len(self.time_flyby_t[int_][i])):
                        if self.freq_flyby_t[int_][i][k] > 0:
                            time_t_arr.append(self.time_flyby_t[int_][i][k])
                            freq_t_arr.append(self.freq_flyby_t[int_][i][k])
                            marker_t_arr.append(self.fb_t_SMBH[int_][i][k])

            marker_nn_arr = np.asarray(marker_nn_arr)
            marker_t_arr = np.asarray(marker_t_arr)

            xdata_set1[int_].append(time_nn_arr)
            xdata_set1[int_] = np.asarray(xdata_set1[int_])
            ydata_set1[int_].append(freq_nn_arr)
            ydata_set1[int_] = np.asarray(ydata_set1[int_])
            xdata_set2[int_].append(time_t_arr)
            xdata_set2[int_] = np.asarray(xdata_set2[int_])
            ydata_set2[int_].append(freq_t_arr)
            ydata_set2[int_] = np.asarray(ydata_set2[int_])
            
            for i in range(2):
                plot_init.tickers(ax[i], 'plot')
                ax[i].set_xlabel('Time [Myr]')
                ax[i].set_ylabel(r' $\log_{10} f_{\rm{GW}}$ [Hz]')
                ax[0].set_xlim(0, 1.1*max(xdata_set1[int_][0]))
                ax[1].set_xlim(0, 1.1*max(xdata_set2[int_][0]))
                if i == 0:
                    ax[i].set_title(plot_title[i])
                    ax[i].scatter(xdata_set1[int_][0][marker_nn_arr > 0], 
                                  np.log10(ydata_set1[int_][0][marker_nn_arr > 0]), 
                                  edgecolors = 'black', color = colours[int_], 
                                  marker = 'X')
                    ax[i].scatter(xdata_set1[int_][0][marker_nn_arr < 1], 
                                  np.log10(ydata_set1[int_][0][marker_nn_arr < 1]), 
                                  edgecolors = 'black', color = colours[int_], 
                                  label = self.integrator[int_])
                else:
                    ax[i].set_title(plot_title[i])
                    ax[i].scatter(xdata_set2[int_][0][marker_t_arr[int_] > 0], 
                                  np.log10(ydata_set2[int_][0][marker_t_arr[int_] > 0]), 
                                  edgecolors = 'black', color = colours[int_], 
                                  marker = 'X')
                    ax[i].scatter(xdata_set2[int_][marker_t_arr[int_] < 1], 
                                  np.log10(ydata_set2[int_][marker_t_arr[int_] < 1]), 
                                  edgecolors = 'black', color = colours[int_])

        ax1.legend()
        plt.savefig('figures/gravitational_waves/events_time.pdf', dpi = 300, bbox_inches='tight')
