from amuse.lab import *
from file_logistics import *
from scipy.special import jv 
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
        self.tH = (self.H0*(3.2408*10**-20))**-1 * 1/(3600*365*24*10**6) | units.Myr
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

                        mass_IMBH = []
                        semi_SMBH_GW_indiv = []
                        semi_NN_GW_indiv = []
                        semi_t_GW_indiv = []
                        ecc_SMBH_GW_indiv = []
                        ecc_NN_GW_indiv = []
                        ecc_t_GW_indiv = []
                        nharm_SMBH_GW_indiv = []
                        nharm_NN_GW_indiv = []
                        nharm_t_GW_indiv = []

                        mass1 = data.iloc[parti_][0][1]
                        SMBH_sys_mass = data.iloc[0][0][1]
                        Nclose_indiv = 0
                        Nfbnn_indiv = 0
                        Nfbt_indiv = 0

                        freq_SMBH_GW_indiv = [-5]
                        freq_NN_GW_indiv = [-5]
                        freq_t_GW_indiv = [-5]
                        strain_SMBH_GW_indiv = [-5]
                        strain_NN_GW_indiv = [-5]
                        strain_t_GW_indiv = [-5]
                        time_SMBH_GW_indiv = [-5]
                        time_NN_GW_indiv = [-5]
                        time_t_GW_indiv = [-5]
                        SMBH_NN_event = [-5]
                        SMBH_t_event = [-5]

                        if parti_ != 0:
                            for col_ in range(np.shape(data)[1]-1):
                                sem_SMBH = data.iloc[parti_][col_][7][0]
                                ecc_SMBH = data.iloc[parti_][col_][8][0]

                                strain_SMBH = self.gw_strain(sem_SMBH, ecc_SMBH, mass1, SMBH_sys_mass)
                                freq_SMBH = self.gw_freq(sem_SMBH, ecc_SMBH, mass1, SMBH_sys_mass)
                                nharm_SMBH = self.gw_harmonic_mode(ecc_SMBH)
                                if strain_SMBH > 10**-27 and freq_SMBH > 10**-12:
                                    strain_SMBH_GW_indiv.append(strain_SMBH)
                                    freq_SMBH_GW_indiv.append(freq_SMBH)
                                    nharm_SMBH_GW_indiv.append(nharm_SMBH)
                                    time_SMBH_GW_indiv.append(10**-3 * col_)
                                    semi_SMBH_GW_indiv.append(sem_SMBH)
                                    ecc_SMBH_GW_indiv.append(ecc_SMBH)
                                
                                semi_major_nn = abs(data.iloc[parti_][col_][7][1])
                                semi_major_t = abs(data.iloc[parti_][col_][7][2])
                                ecc_nn = (data.iloc[parti_][col_][8][1])
                                ecc_t = (data.iloc[parti_][col_][8][2])
                                for part_ in range(np.shape(data)[0]):
                                    if data.iloc[part_][0][0] == data.iloc[parti_][col_][6][1]:
                                        mass2 = data.iloc[part_][0][1]

                                strain_nn = self.gw_strain(semi_major_nn, ecc_nn, mass1, mass2)
                                freq_nn = self.gw_freq(semi_major_nn, ecc_nn, mass1, mass2)
                                nharm_nn = self.gw_harmonic_mode(ecc_nn)
                                if strain_nn > 10**-27 and freq_nn > 10**-12:
                                    strain_NN_GW_indiv.append(strain_nn)
                                    freq_NN_GW_indiv.append(freq_nn)
                                    nharm_NN_GW_indiv.append(nharm_nn)
                                    time_NN_GW_indiv.append(10**-3 * col_)
                                    semi_NN_GW_indiv.append(semi_major_nn)
                                    ecc_NN_GW_indiv.append(ecc_nn)

                                    linex = data.iloc[parti_][col_][2][0] - data.iloc[0][col_][2][0]
                                    liney = data.iloc[parti_][col_][2][1] - data.iloc[0][col_][2][1]
                                    linez = data.iloc[parti_][col_][2][2] - data.iloc[0][col_][2][2]
                                    dist_SMBH = (linex**2+liney**2+linez**2).sqrt()
                                    dist_NN = data.iloc[parti_][col_][-1]

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

                                strain_t = self.gw_strain(semi_major_t, ecc_t, mass1, mass2)
                                freq_t = self.gw_freq(semi_major_t, ecc_t, mass1, mass2)
                                nharm_t = self.gw_harmonic_mode(ecc_t)
                                if strain_t > 10**-27 and freq_t > 10**-12:
                                    strain_t_GW_indiv.append(strain_t)
                                    freq_t_GW_indiv.append(freq_t)
                                    nharm_t_GW_indiv.append(nharm_t)
                                    time_t_GW_indiv.append(10**-3 * col_) 
                                    semi_t_GW_indiv.append(semi_major_t)
                                    ecc_t_GW_indiv.append(ecc_t)    

                                    if (tSMBH):
                                        Nfbt_indiv += 1
                                        SMBH_t_event.append(1)
                                    else:
                                        Nfbt_indiv += 0.5
                                        SMBH_t_event.append(-5)

                                if data.iloc[parti_][col_][-1] < 5e-3:
                                    Nclose_indiv += 1

                            mass_IMBH = np.asarray(mass_IMBH)
                            Ntot_indiv = Nfbnn_indiv + Nfbt_indiv

                            path = '/home/erwanh/Desktop/SteadyStateBH/Data_Process/data/tGW/'
                            stab_tracker = pd.DataFrame()
                            df_stabtime = pd.Series({'Integrator': self.integrator[int_],
                                                     'Simulation Time': 10**3*sim_time,
                                                     'Population': 10*round(0.1*np.shape(data)[0]),
                                                     'mass SMBH': SMBH_sys_mass,
                                                     'mass IMBH1': mass1,
                                                     'mass IMBH2': mass2,
                                                     'No. Binary Events': Nfbnn_indiv,
                                                     'No. Tertiary Events': Nfbt_indiv,
                                                     'No. Total Events': Ntot_indiv,
                                                     'No. Close Encounter': Nclose_indiv,
                                                     'SMBH Binary Frequencies': freq_SMBH_GW_indiv,
                                                     'SMBH Binary Strain': strain_SMBH_GW_indiv,
                                                     'SMBH Binary Time': time_SMBH_GW_indiv,
                                                     'SMBH Binary Semi major': semi_SMBH_GW_indiv,
                                                     'SMBH Binary Eccentricity': ecc_SMBH_GW_indiv,
                                                     'FlyBy Binary Frequencies': freq_NN_GW_indiv,
                                                     'Flyby Binary Strain': strain_NN_GW_indiv,
                                                     'Flyby Binary Time': time_NN_GW_indiv,
                                                     'Flyby Binary Semi major': semi_NN_GW_indiv,
                                                     'Flyby Binary Eccentricity': ecc_NN_GW_indiv,
                                                     'Flyby SMBH Event': SMBH_NN_event,
                                                     'Flyby Tertiary Frequencies': freq_t_GW_indiv,
                                                     'Flyby Tertiary Strain': strain_t_GW_indiv,
                                                     'Flyby Tertiary Time': time_t_GW_indiv,
                                                     'Flyby Tertiary Semi major': semi_t_GW_indiv,
                                                     'Flyby Tertiary Eccentricity': ecc_t_GW_indiv,
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
        self.mass_parti = [[ ], [ ]]
        self.mass_IMBH = [[ ], [ ]]

        self.freq_flyby_SMBH = [[ ], [ ]]
        self.strain_flyby_SMBH = [[ ], [ ]]
        self.time_flyby_SMBH = [[ ], [ ]]
        self.semi_flyby_SMBH = [[ ], [ ]]
        self.ecc_flyby_SMBH = [[ ], [ ]]

        self.freq_flyby_nn = [[ ], [ ]]
        self.strain_flyby_nn = [[ ], [ ]]
        self.time_flyby_nn = [[ ], [ ]]
        self.semi_flyby_nn = [[ ], [ ]]
        self.ecc_flyby_nn = [[ ], [ ]]

        self.freq_flyby_t = [[ ], [ ]]
        self.strain_flyby_t = [[ ], [ ]]
        self.time_flyby_t = [[ ], [ ]]
        self.semi_flyby_t = [[ ], [ ]]
        self.ecc_flyby_t = [[ ], [ ]]

        self.tot_sim_time = [[ ], [ ]]
        self.tot_events = [[ ], [ ]]
        self.fb_nn_events = [[ ], [ ]]
        self.fb_nn_SMBH = [[ ], [ ]]
        self.fb_t_events = [[ ], [ ]]
        self.fb_t_SMBH = [[ ], [ ]]

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
                self.mass_SMBH[int_idx].append(data_file.iloc[0][3])
                self.mass_parti[int_idx].append(data_file.iloc[0][4])
                self.mass_IMBH[int_idx].append(data_file.iloc[0][5])
                self.fb_nn_events[int_idx].append(data_file.iloc[0][6])
                self.fb_t_events[int_idx].append(data_file.iloc[0][7])
                self.tot_events[int_idx].append(data_file.iloc[0][8])
                self.close_enc[int_idx].append(data_file.iloc[0][9])

                self.freq_flyby_SMBH[int_idx].append(data_file.iloc[0][10])
                self.strain_flyby_SMBH[int_idx].append(data_file.iloc[0][11])
                self.time_flyby_SMBH[int_idx].append(data_file.iloc[0][12])
                self.semi_flyby_SMBH[int_idx].append(data_file.iloc[0][13])
                self.ecc_flyby_SMBH[int_idx].append(data_file.iloc[0][14])

                self.freq_flyby_nn[int_idx].append(data_file.iloc[0][15])
                self.strain_flyby_nn[int_idx].append(data_file.iloc[0][16])
                self.time_flyby_nn[int_idx].append(data_file.iloc[0][17])
                self.semi_flyby_nn[int_idx].append(data_file.iloc[0][18])
                self.ecc_flyby_nn[int_idx].append(data_file.iloc[0][19])
                self.fb_nn_SMBH[int_idx].append(data_file.iloc[0][20])

                self.freq_flyby_t[int_idx].append(data_file.iloc[0][21])
                self.strain_flyby_t[int_idx].append(data_file.iloc[0][22])
                self.time_flyby_t[int_idx].append(data_file.iloc[0][23])
                self.semi_flyby_t[int_idx].append(data_file.iloc[0][24])
                self.ecc_flyby_t[int_idx].append(data_file.iloc[0][25])
                self.fb_t_SMBH[int_idx].append(data_file.iloc[0][26])

        for int_ in range(2):
            self.sim_time[int_] = np.asarray(self.sim_time[int_], dtype = 'object')
            self.pop[int_] = np.asarray(self.pop[int_], dtype = 'object')
            self.mass_SMBH[int_] = np.asarray(self.mass_SMBH[int_], dtype = 'object')
            self.mass_parti[int_] = np.asarray(self.mass_parti[int_], dtype = 'object')
            self.mass_IMBH[int_] = np.asarray(self.mass_IMBH[int_], dtype = 'object')
            self.fb_nn_events[int_] = np.asarray(self.fb_nn_events[int_], dtype = 'object')
            self.fb_t_events[int_] = np.asarray(self.fb_t_events[int_], dtype = 'object')
            self.tot_events[int_] = np.asarray(self.tot_events[int_], dtype = 'object')
            self.close_enc[int_] = np.asarray(self.close_enc[int_], dtype = 'object')
            self.freq_flyby_SMBH[int_] = np.asarray(self.freq_flyby_SMBH[int_], dtype = 'object')
            self.strain_flyby_SMBH[int_] = np.asarray(self.strain_flyby_SMBH[int_], dtype = 'object')
            self.time_flyby_SMBH[int_] = np.asarray(self.time_flyby_SMBH[int_], dtype = 'object')
            self.semi_flyby_SMBH[int_] = np.asarray(self.semi_flyby_SMBH[int_], dtype = 'object')
            self.ecc_flyby_SMBH[int_] = np.asarray(self.ecc_flyby_SMBH[int_], dtype = 'object')
            self.freq_flyby_nn[int_] = np.asarray(self.freq_flyby_nn[int_], dtype = 'object')
            self.strain_flyby_nn[int_] = np.asarray(self.strain_flyby_nn[int_], dtype = 'object')
            self.time_flyby_nn[int_] = np.asarray(self.time_flyby_nn[int_], dtype = 'object')
            self.semi_flyby_nn[int_] = np.asarray(self.semi_flyby_nn[int_], dtype = 'object')
            self.ecc_flyby_nn[int_] = np.asarray(self.ecc_flyby_nn[int_], dtype = 'object')
            self.fb_nn_SMBH[int_] = np.asarray(self.fb_nn_SMBH[int_], dtype = 'object')
            self.freq_flyby_t[int_] = np.asarray(self.freq_flyby_t[int_], dtype = 'object')
            self.strain_flyby_t[int_] = np.asarray(self.strain_flyby_t[int_], dtype = 'object')
            self.time_flyby_t[int_] = np.asarray(self.time_flyby_t[int_], dtype = 'object')
            self.semi_flyby_t[int_] = np.asarray(self.semi_flyby_t[int_], dtype = 'object')
            self.ecc_flyby_t[int_] = np.asarray(self.ecc_flyby_t[int_], dtype = 'object')
            self.fb_t_SMBH[int_] = np.asarray(self.fb_t_SMBH[int_], dtype = 'object')
        
        with open('figures/gravitational_waves/output/event_rate.txt', 'w') as file:
                for int_ in range(1):
                    pop_idx = np.where(self.pop[int_] == 10)
                    tot_events_t = self.tot_events[int_][pop_idx]
                    tot_events = tot_events_t[tot_events_t > 0]
                    
                    obs_event_nn = 0
                    samples_nn = 0
                    perc_SMBH_nn = 0
                    obs_event_t = 0
                    samples_t = 0
                    perc_SMBH_t = 0
                    obs_event_SMBH = 0
                    samples_SMBH = 0
                    no_strain_nn = 0
                    no_strain_t = 0
                    no_strain_S = 0
                    
                    no_data = np.shape(self.freq_flyby_nn[int_][pop_idx])[0]
                    for parti_ in range(no_data):
                        sim_time = self.sim_time[int_][pop_idx][parti_]

                        freq_temp = np.asarray(self.freq_flyby_nn[int_][pop_idx][parti_])
                        strain_temp = np.asarray(self.strain_flyby_nn[int_][pop_idx][parti_])
                        fbb_IMBH = np.asarray(self.fb_nn_SMBH[int_][pop_idx][parti_])
                        obs_fbnn = freq_temp[(freq_temp > 10**-4) & (fbb_IMBH > 0) & (strain_temp > 10**-22)]
                        obs_event_nn += (len(obs_fbnn)/sim_time)
                        perc_SMBH_nn += len(obs_fbnn)/len(freq_temp)
                        if len(obs_fbnn) > 0:
                            samples_nn += 1
                            no_strain_nn += len(obs_fbnn) - len(freq_temp[(freq_temp > 10**-4) & (fbb_IMBH > 0)])
                        
                        freq_temp = np.asarray(self.freq_flyby_t[int_][pop_idx][parti_])
                        strain_temp = np.asarray(self.strain_flyby_t[int_][pop_idx][parti_])
                        fbt_IMBH = np.asarray(self.fb_t_SMBH[int_][pop_idx][parti_])
                        obs_fbt = freq_temp[(freq_temp > 10**-4) & (fbt_IMBH > 0) & (strain_temp > 10**-22)]
                        obs_event_t += (len(obs_fbt)/sim_time)
                        perc_SMBH_t += len(obs_fbt)/len(freq_temp)
                        if len(obs_fbt) > 0:
                            samples_t += 1
                            no_strain_t += len(obs_fbt) - len(freq_temp[(freq_temp > 10**-4) & (fbt_IMBH > 0)])

                        freq_temp = np.asarray(self.freq_flyby_SMBH[int_][pop_idx][parti_])
                        strain_temp = np.asarray(self.strain_flyby_SMBH[int_][pop_idx][parti_])
                        obs_fbS = freq_temp[(freq_temp > 10**-4) & (strain_temp > 10**-22)]
                        obs_event_SMBH += (len(obs_fbS)-len(obs_fbt)-len(obs_fbnn))/sim_time
                        if len(obs_fbS) > 0:
                            samples_SMBH += 1
                            no_strain_S += len(obs_fbS) - len(freq_temp[(freq_temp > 10**-4)])

                    obs_event_nn /= max(1, samples_nn)
                    obs_event_t /= max(1, samples_t)
                    obs_event_SMBH /= max(1, samples_SMBH)
                    no_strain_nn /= no_data
                    no_strain_t /= no_data
                    no_strain_S /= no_data
                    perc_SMBH_nn = (obs_event_SMBH)/(obs_event_SMBH+obs_event_nn)
                    perc_SMBH_t = (obs_event_t)/(obs_event_SMBH+obs_event_t)

                    fb_nn_events = self.fb_nn_events[int_]#[pop_idx]
                    fb_nn_events = fb_nn_events[tot_events_t > 0]
                    fb_t_events = self.fb_t_events[int_]#[pop_idx]
                    fb_t_events = fb_t_events[tot_events_t > 0]
                    sim_time = self.sim_time[int_]#[pop_idx]
                    sim_time = sim_time[tot_events_t > 0]

                    if len(tot_events) > 0:
                        if int_ == 0:
                            file.write('Detectable means f > 1e-4 Hz and strain > 1e-22 when observing a source 1Mpc away')
                        file.write('\nData for '+str(self.integrator[int_]))
                        file.write('\nAverage total event rate per Myr:                      ' + str(np.mean(tot_events/sim_time * 10**6)))
                        file.write('\nAverage total flyby events per Myr:                    ' + str(np.mean(fb_nn_events/sim_time) * 10**6))
                        file.write('\nAverage total tertiary events per Myr:                 ' + str(np.mean(fb_t_events/sim_time) * 10**6))
                        file.write('\nAverage detectable IMBH-IMBH per Myr:                  ' + str((obs_event_nn + obs_event_t) * 10**6))
                        file.write('\nAverage detectable flyby IMBH-IMBH per Myr:            ' + str(obs_event_nn * 10**6))
                        file.write('\nAverage detectable tertiary IMBH-IMBH per Myr:         ' + str(obs_event_t * 10**6))
                        file.write('\nAverage detectable SMBH-IMBH event rate per Myr:       ' + str(obs_event_SMBH * 10**6))
                        file.write('\nPercentage SMBH flyby events:                          ' + str(perc_SMBH_nn))
                        file.write('\nPercentage SMBH tertiary events:                       ' + str(perc_SMBH_t))
                        file.write('\nAvg. fraction in freq but not detectable IMBH-IMBH:    ' + str((no_strain_nn+no_strain_t)/2))
                        file.write('\nAvg. fraction in freq but not detectable SMBH-IMBH:    ' + str(no_strain_S))
                        file.write('\n========================================================================')

    def coll_radius(self, mass_arr):
        return 3 * (2*constants.G*mass_arr)/(constants.c**2)

    def forecast_interferometer(self, ax, m1, m2):
        """
        Function to plot the LISA and aLIGO frequency range in Ge a vs. (1-e) plots
        """
        
        ecc_range = np.linspace(0.0001, (1-10**-8), 50)

        self.LIGO_semimaj_max = self.gw_cfreq_semi(ecc_range[1:], 10000 | units.Hz, m1, m2)
        self.LIGO_semimaj_min = self.gw_cfreq_semi(ecc_range[1:], 10 | units.Hz, m1, m2)
        self.LIGO_semimaj = self.gw_cfreq_semi(ecc_range[1:], 200 | units.Hz, m1, m2)

        self.LISA_semimaj_max = self.gw_cfreq_semi(ecc_range[1:], 1 | units.Hz, m1, m2)
        self.LISA_semimaj_min = self.gw_cfreq_semi(ecc_range[1:], 1e-4 | units.Hz, m1, m2)
        self.LISA_semimaj = self.gw_cfreq_semi(ecc_range[1:], 1e-2 | units.Hz, m1, m2)

        ecc_range = [np.log(1-i) for i in ecc_range[1:]]
        self.text_angle = np.degrees(np.arctan((ecc_range[30]-ecc_range[20])/(self.LIGO_semimaj[30]-self.LIGO_semimaj[20]))) - 20

        ax.plot(self.LIGO_semimaj_min, ecc_range, linestyle = ':', color = 'white')
        ax.plot(self.LIGO_semimaj, ecc_range, linestyle = '-.', color = 'white')
        ax.plot(self.LIGO_semimaj_max, ecc_range, linestyle = ':', color = 'white')
        ax.fill_between(np.append(self.LIGO_semimaj_min, self.LIGO_semimaj_max[::-1]), 
                        np.append(ecc_range[:], ecc_range[::-1]), alpha = 0.6, color = 'blue')
        ax.plot(self.LISA_semimaj_min, ecc_range, linestyle = ':', color = 'white')
        ax.plot(self.LISA_semimaj, ecc_range, linestyle = '-.', color = 'white')
        ax.plot(self.LISA_semimaj_max, ecc_range, linestyle = ':', color = 'white')
        ax.fill_between(np.append(self.LISA_semimaj_min, self.LISA_semimaj_max[::-1]), 
                        np.append(ecc_range[:], ecc_range[::-1]), alpha = 0.6, color = 'red')

        return ax

    def gfunc(self, ecc):
        nharm = self.gw_harmonic_mode(ecc)
        return nharm**4/32 * ((jv(nharm-2, nharm*ecc)-2*ecc*jv(nharm-1, nharm*ecc) + 2/nharm * jv(nharm, nharm*ecc) + 2*ecc*jv(nharm+1, nharm*ecc) - jv(nharm+2, nharm*ecc))**2 + (1-ecc**2)*(jv(nharm-2, nharm*ecc) - 2*jv(nharm, nharm*ecc) + jv(nharm+2, nharm*ecc))**2 + 4/(3*nharm**2)*(jv(nharm, nharm*ecc)**2))
    
    def gw_cfreq_semi(self, ecc_arr, freq_val, m1, m2):
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

    def gw_freq(self, semi, ecc, m1, m2):
        """
        Frequency equation is based on Samsing et al. 2014 eqn (43). 
        
        Inputs:
        semi:   The semi-major axes of the system
        ecc:    The eccentricity of the binary system
        m1/m2:  The individual mass components of the binary system
        """

        nharm = self.gw_harmonic_mode(ecc)
        freq =  (2*np.pi)**-1*np.sqrt(constants.G*(m1+m2)/abs(semi)**3) * nharm
        return freq.value_in(units.Hz)

    def gw_harmonic_mode(self, ecc):
        """
        Finding the peak harmonic of gravitational frequency for a given eccentric orbit.
        Equation 36 of Wen (2003)
        
        Inputs:
        ecc:    The eccentricity of the orbit
        """
        nharm = 2*(1+ecc)**1.1954/(1-ecc**2)**1.5
        return nharm

    def gw_strain(self, semi, ecc, m1, m2):
        """
        Use of eqn (7) Kremer et al. 2018.
        Use of eqn (20) of Peters and Matthews (1963).
        At 1Mpc, z ~ 0.0228 from cosmology calc.
        
        Inputs:
        semi:   The semi-major axes of the system
        ecc:    The eccentricity of the binary system
        m1/m2:  The individual mass components of the binary system
        """

        dist = 1 | units.Mpc  # Taken from [https://imagine.gsfc.nasa.gov/features/cosmic/milkyway_info.html]
        redshift = 0.0228

        chirp_mass = (m1*m2)**0.6/(m1+m2)**0.2 * (1+redshift)**-1
        cfactor = 2/(3*np.pi**(4/3)) * (constants.G**(5/3))/(constants.c**3) * (dist*(1+redshift))**-2
        ecc_func = (1+(73/24)*ecc**2+(37/96)*ecc**4)*(1-ecc**2)**-3.5

        nharm = self.gw_harmonic_mode(ecc)
        freq = self.gw_freq(semi, ecc, m1, m2) * (1 | units.Hz)

        strain = cfactor * chirp_mass**(5/3) * freq**(-1/3) * (2/nharm)**(2/3) * (self.gfunc(ecc)/ecc_func)
        return (strain.value_in(units.s**-1.6653345369377348e-16))**0.5

    def gw_timescale(self, semi, ecc, m1, m2):
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
        ax.fill_between(PTA_freq, -30, -13, alpha = 0.35, color = 'red')
        ax.text(-2.3, -20.5, 'Lisa Sensitivity', fontsize ='small', rotation = 276)
        ax.text(-6, -15, 'PTA Window', fontsize ='small', rotation = 270)
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

        tgw_nn = [[ ], [ ]]
        semi_nn = [[ ], [ ]]
        ecc_nn = [[ ], [ ]]
        freqnn_det = [[ ], [ ]]

        tgw_t = [[ ], [ ]]
        semi_t = [[ ], [ ]]
        ecc_t = [[ ], [ ]]
        freqt_det = [[ ], [ ]]
        
        for int_ in range(1):
            for parti_ in range(len(self.semi_flyby_nn[int_])): #Looping through every individual particle
                for event_ in range(len(self.semi_flyby_nn[int_][parti_])): #Looping through every detected event
                    semi_fb_nn = self.semi_flyby_nn[int_][parti_][event_]
                    if np.asarray(self.fb_nn_SMBH[int_][parti_][event_]) < 0 and semi_fb_nn < 1 | units.parsec:
                        ecc_fb_nn = self.ecc_flyby_nn[int_][parti_][event_]
                        gwtime_nn = self.gw_timescale(semi_fb_nn, ecc_fb_nn, self.mass_IMBH[int_][parti_], 
                                                      self.mass_parti[int_][parti_]).value_in(units.Myr)

                        tgw_nn[int_].append(gwtime_nn)
                        semi_nn[int_].append(semi_fb_nn.value_in(units.pc))
                        ecc_nn[int_].append(np.log10(1-ecc_fb_nn))

            for parti_ in range(len(self.semi_flyby_t[int_])):
                for event_ in range(len(self.semi_flyby_t[int_][parti_])):
                    semi_fb_t = self.semi_flyby_t[int_][parti_][event_]
                    if np.asarray(self.fb_t_SMBH[int_][parti_][event_]) < 0 and semi_fb_t < 1 | units.parsec:
                        ecc_fb_t = self.ecc_flyby_t[int_][parti_][event_]
                        gwtime_t = self.gw_timescale(semi_fb_t, ecc_fb_t, self.mass_parti[int_][parti_], 
                                                     self.mass_IMBH[int_][parti_]).value_in(units.Myr)

                        tgw_t[int_].append(gwtime_t)
                        semi_t[int_].append(semi_fb_t.value_in(units.pc))
                        ecc_t[int_].append(np.log10(1-ecc_fb_t))

                        if self.freq_flyby_t[int_][parti_][event_] > 10**-11:
                            freqt_det[int_].append(1)
                        else:
                            freqt_det[int_].append(0)

            tgw_nn[int_] = np.asarray(tgw_nn[int_])
            ecc_nn[int_] = np.asarray(ecc_nn[int_])
            semi_nn[int_] = np.asarray(semi_nn[int_])  
            tgw_t[int_] = np.asarray(tgw_t[int_])
            ecc_t[int_] = np.asarray(ecc_t[int_])
            semi_t[int_] = np.asarray(semi_t[int_])
        
        ############### PLOTTING OF a vs. (1-e) FOR BIN ##############

        xmin = -10
        xmax = 0
        ymin = -4

        red_mass = (self.mass_parti[0][0]**2)*(2*self.mass_parti[0][0])
        x_arr = np.linspace(-10, 0, 2000)
        const_tgw = [np.log10(1-np.sqrt(1-((256*self.tH*(constants.G**3)/(5*constants.c**5)*red_mass*(10**(i) * (1 | units.pc)) **-4)) **(1/3.5))) for i in x_arr]

        fig = plt.figure(figsize=(15, 6))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax = [ax1, ax2]
        extent = [1.1*xmin, 1.1*xmax, 1.1*ymin, 0]
        for j in range(1):
            bin2d_sim, xedg, xedg, image = ax[j].hist2d(np.log10(semi_nn[j]), ecc_nn[j], bins=(200, 200), 
                                                      range=([1.1*xmin, 1.1*xmax], [1.1*ymin, 0]))
            bin2d_sim, xedg, xedg, image = ax[j].hist2d(np.log10(semi_t[j]), ecc_t[j], bins=(200, 200), 
                                                      range=([1.1*xmin, 1.1*xmax], [1.1*ymin, 0]))
                                                      
            bin2d_sim /= np.max(bin2d_sim)
            contours = ax[j].imshow((bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
            
            ax[j].set_ylabel(r'$\log_{10}(1-e)$')
            ax[j].set_xlabel(r'$\log_{10} a$ [pc]')
            ax[j].set_xlim(xmin, xmax)
            ax[j].set_ylim(ymin, 0)
                
            plot_init.tickers(ax[j], 'hist')
            """if j == 0:
                ax1.scatter((semi_nn[j][tgw_nn[j] > self.tH]), ecc_nn[j][tgw_nn[j] > self.tH], s = 8, color = 'black', label = 'Nearest Neighbour')
                ax1.scatter((semi_t[j][tgw_t[j] > self.tH]), ecc_t[j][tgw_t[j] > self.tH], s = 8, color = 'purple', label = 'Tertiary')
                ax1.scatter((semi_nn[j][tgw_nn[j] < self.tH]), ecc_nn[j][tgw_nn[j] < self.tH], s = 8, color = 'black', marker = 'X', label = r'$t_{\rm{GW}} < t_H$')
                ax1.scatter((semi_t[j][tgw_t[j] < self.tH]), ecc_t[j][tgw_t[j] < self.tH], s = 8, color = 'purple', marker = 'X')

            else:
                ax1.scatter((semi_nn[j][tgw_nn[j] > self.tH]), ecc_nn[j][tgw_nn[j] > self.tH], s = 8, color = 'black')
                ax1.scatter((semi_t[j][tgw_t[j] > self.tH]), ecc_t[j][tgw_t[j] > self.tH], s = 8, color = 'purple')
                ax1.scatter((semi_nn[j][tgw_nn[j] < self.tH]), ecc_nn[j][tgw_nn[j] < self.tH], s = 8, color = 'black', marker = 'X')
                ax1.scatter((semi_t[j][tgw_t[j] < self.tH]), ecc_t[j][tgw_t[j] < self.tH], s = 8, color = 'purple', marker = 'X')"""

            self.forecast_interferometer(ax[j], self.mass_parti[j][0], self.mass_IMBH[j][0])
            ax[j].plot(x_arr, const_tgw, color = 'white')
            ax[j].text(-9.5, -3, r'aLIGO ($f_{\rm{peak}} = 200$ Hz)', verticalalignment = 'center', fontsize ='small', rotation=self.text_angle+8.7, color = 'white')
            ax[j].text(-6.6, -3, r'LISA ($f_{\rm{peak}} = 10^{-2}$ Hz)', verticalalignment = 'center', fontsize ='small', rotation=self.text_angle+8.7, color = 'white')
            ax[j].text(-2.55, -3, r'$t_{\rm{GW}} > t_H$', verticalalignment = 'center', fontsize ='small', rotation=self.text_angle+18.5, color = 'white')
            ax[j].text(-2.95, -3, r'$t_{\rm{GW}} < t_H$', verticalalignment = 'center', fontsize ='small', rotation=self.text_angle+18.5, color = 'white')

        plt.savefig('figures/gravitational_waves/ecc_semi_bins_IMBH_histogram.pdf', dpi=300, bbox_inches='tight')
        plt.clf()

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
                grav_avg_time = self.gw_timescale(self.semi_SMBH_avg[int_][i] * (1 | units.pc), self.ecc_SMBH_avg[int_][i], 
                                             self.mass_SMBH[int_][i][0], self.mass_SMBH[int_][i][1]).value_in(units.Myr)
                tgw_SMBH_avg[int_].append(grav_avg_time)
                grav_min_time = self.gw_timescale(self.semi_SMBH_min[int_][i] * (1 | units.pc), self.ecc_SMBH_min[int_][i], 
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
            self.forecast_interferometer(ax1, self.mass_SMBH[int_][0], self.mass_parti[int_][0])
            plt.colorbar(colour_axes, ax = ax1, label = r'$\log_{10} \langle t_{GW}\rangle$ [Myr]')
            ax1.set_xlim(-10, 1.1*max(xmax_ecc_sem))
            ax1.text(-8.95, -2, r'aLIGO ($f_{\rm{peak}} = 200$ Hz)', verticalalignment = 'center', fontsize ='small', rotation=self.text_angle+9.05, color = 'white')
            ax1.text(-6.05, -2, r'LISA ($f_{\rm{peak}} = 10^{-2}$ Hz)', verticalalignment = 'center', fontsize ='small', rotation=self.text_angle+9.05, color = 'white')
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
            ax.set_ylim(-30, -13.2)
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

print('...tGW_plotters...')
cst = coupled_systems()
#cst.new_data_extractor()
cst.combine_data()
cst.IMBH_tgw_plotter()
cst.SMBH_tgw_plotter()
cst.strain_freq_plotter()
cst.transient_events()