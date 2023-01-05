from amuse.lab import *
from file_logistics import *
from scipy.special import jv 
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import warnings
import pandas as pd
import statsmodels.api as sm
import LISA_Curves.LISA as li
import pickle as pkl

class gw_calcs(object):
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
        Hermite_data = glob.glob(os.path.join('/media/erwanh/Elements/Hermite/particle_trajectory_tGW_new/*'))
        chaoticH = ['data/Hermite/no_addition/chaotic_simulation/'+str(i[59:]) for i in Hermite_data]
        GRX_data = glob.glob(os.path.join('/media/erwanh/Elements/GRX/particle_trajectory_tGW_new/*'))
        chaoticG = ['data/GRX/no_addition/chaotic_simulation/'+str(i[55:]) for i in GRX_data]

        filename = [natsort.natsorted(Hermite_data), natsort.natsorted(GRX_data)] 
        filenameC = [natsort.natsorted(chaoticH), natsort.natsorted(chaoticG)] #CRASHED AT SIM30, 181
        
        for int_ in range(2):
            for file_ in range(len(filename[int_])):
                with open(filenameC[int_][file_], 'rb') as input_file:
                    chaotic_tracker = pkl.load(input_file)
                    if chaotic_tracker.iloc[0][6] <= 40 and chaotic_tracker.iloc[0][6] > 5:
                        with open(filename[int_][file_], 'rb') as input_file:
                            print('Reading file', file_, ': ', filename[int_][file_])
                            data = pkl.load(input_file)
                            sim_time = np.shape(data)[1]-1
                            SMBH_sys_mass = data.iloc[0][0][1]

                            for parti_ in range(np.shape(data)[0]):
                                count = len(fnmatch.filter(os.listdir('data/tGW/'), '*.*'))
                                mass1 = data.iloc[parti_][0][1]
                                tot_mass = SMBH_sys_mass + mass1
                                if not (isinstance(data.iloc[parti_][-1][0], np.uint64)) or data.iloc[parti_][-1][1] == tot_mass:
                                    merge_Bool = 1
                                else: 
                                    merge_Bool = 0

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
                                        sem_SMBH = abs(data.iloc[parti_][col_][7][0])
                                        ecc_SMBH = data.iloc[parti_][col_][8][0]

                                        strain_SMBH = self.gw_strain(sem_SMBH, ecc_SMBH, mass1, SMBH_sys_mass)
                                        freq_SMBH = self.gw_freq(sem_SMBH, ecc_SMBH, mass1, SMBH_sys_mass)
                                        nharm_SMBH = self.gw_harmonic_mode(ecc_SMBH)
                                        if freq_SMBH > 10**-12:
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
                                        if freq_nn > 10**-12:
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

                                            if dist_SMBH.value_in(units.pc) == dist_NN or mass2 > 10**5 | units.MSun:
                                                Nfbnn_indiv += 1
                                                SMBH_NN_event.append(1)
                                            else:
                                                Nfbnn_indiv += 0.5
                                                SMBH_NN_event.append(-5)
                                                
                                        tSMBH = False
                                        if sem_SMBH == semi_major_t or ecc_SMBH == ecc_t:
                                            mass2 = SMBH_sys_mass
                                            tSMBH = True

                                        strain_t = self.gw_strain(semi_major_t, ecc_t, mass1, mass2)
                                        freq_t = self.gw_freq(semi_major_t, ecc_t, mass1, mass2)
                                        nharm_t = self.gw_harmonic_mode(ecc_t)
                                        if freq_t > 10**-12:
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
                                                            'Tertiary SMBH Event': SMBH_t_event,
                                                            'Merger Boolean': merge_Bool})
                                    stab_tracker = stab_tracker.append(df_stabtime, ignore_index = True)
                                    stab_tracker.to_pickle(os.path.join(path, 'IMBH_'+str(self.integrator[int_])+'_tGW_data_indiv_parti_'+str(count)+'_'+str(parti_)+'_local2.pkl'))

    def combine_data(self, filter, integ):
        """
        Function which extracts ALL data and provides 
        the merger rate of the system.
        """

        self.sim_time = [ ]
        self.pop = [ ]
        self.close_enc = [ ]
        self.mass_SMBH = [ ]
        self.mass_parti = [ ]
        self.mass_IMBH = [ ]

        self.freq_flyby_SMBH = [ ]
        self.strain_flyby_SMBH = [ ]
        self.time_flyby_SMBH = [ ]
        self.semi_flyby_SMBH = [ ]
        self.ecc_flyby_SMBH = [ ]

        self.freq_flyby_nn = [ ]
        self.strain_flyby_nn = [ ]
        self.time_flyby_nn = [ ]
        self.semi_flyby_nn = [ ]
        self.ecc_flyby_nn = [ ]

        self.freq_flyby_t = [ ]
        self.strain_flyby_t = [ ]
        self.time_flyby_t = [ ]
        self.semi_flyby_t = [ ]
        self.ecc_flyby_t = [ ]

        self.tot_sim_time = [ ]
        self.tot_events = [ ]
        self.fb_nn_events = [ ]
        self.fb_nn_SMBH = [ ]
        self.fb_t_events = [ ]
        self.fb_t_SMBH = [ ]

        if integ == 'Hermite':
            int_idx = 0
        else:
            int_idx = 1

        print('Extracting data for ', self.integrator[int_idx])
        files = 0
        tGW_data = natsort.natsorted(glob.glob('data/tGW/*'))
        if filter == 'pop_filt':
            for file_ in range(len(tGW_data)):
                with open(tGW_data[file_], 'rb') as input_file:
                    data_file = pkl.load(input_file)
                    if data_file.iloc[0][2] == 15:
                        if data_file.iloc[0][0] == self.integrator[int_idx]:
                            files += 1
                            self.sim_time.append(data_file.iloc[0][1])
                            self.pop.append(int(data_file.iloc[0][2]))
                            self.mass_SMBH.append(data_file.iloc[0][3])
                            self.mass_parti.append(data_file.iloc[0][4])
                            self.mass_IMBH.append(data_file.iloc[0][5])
                            self.fb_nn_events.append(data_file.iloc[0][6])
                            self.fb_t_events.append(data_file.iloc[0][7])
                            self.tot_events.append(data_file.iloc[0][8])
                            self.close_enc.append(data_file.iloc[0][9])

                            self.freq_flyby_SMBH.append(data_file.iloc[0][10])
                            self.strain_flyby_SMBH.append(data_file.iloc[0][11])
                            self.time_flyby_SMBH.append(data_file.iloc[0][12])
                            self.semi_flyby_SMBH.append(data_file.iloc[0][13])
                            self.ecc_flyby_SMBH.append(data_file.iloc[0][14])

                            self.freq_flyby_nn.append(data_file.iloc[0][15])
                            self.strain_flyby_nn.append(data_file.iloc[0][16])
                            self.time_flyby_nn.append(data_file.iloc[0][17])
                            self.semi_flyby_nn.append(data_file.iloc[0][18])
                            self.ecc_flyby_nn.append(data_file.iloc[0][19])
                            self.fb_nn_SMBH.append(data_file.iloc[0][20])

                            self.freq_flyby_t.append(data_file.iloc[0][21])
                            self.strain_flyby_t.append(data_file.iloc[0][22])
                            self.time_flyby_t.append(data_file.iloc[0][23])
                            self.semi_flyby_t.append(data_file.iloc[0][24])
                            self.ecc_flyby_t.append(data_file.iloc[0][25])
                            self.fb_t_SMBH.append(data_file.iloc[0][26])

        if filter == 'integrator':
            for file_ in range(len(tGW_data)):
                with open(tGW_data[file_], 'rb') as input_file:
                    data_file = pkl.load(input_file)
                    if data_file.iloc[0][0] == 'GRX':
                        int_idx = 1

                        self.sim_time.append(data_file.iloc[0][1])
                        self.pop.append(int(data_file.iloc[0][2]))
                        self.mass_SMBH.append(data_file.iloc[0][3])
                        self.mass_parti.append(data_file.iloc[0][4])
                        self.mass_IMBH.append(data_file.iloc[0][5])
                        self.fb_nn_events.append(data_file.iloc[0][6])
                        self.fb_t_events.append(data_file.iloc[0][7])
                        self.tot_events.append(data_file.iloc[0][8])
                        self.close_enc.append(data_file.iloc[0][9])

                        self.freq_flyby_SMBH.append(data_file.iloc[0][10])
                        self.strain_flyby_SMBH.append(data_file.iloc[0][11])
                        self.time_flyby_SMBH.append(data_file.iloc[0][12])
                        self.semi_flyby_SMBH.append(data_file.iloc[0][13])
                        self.ecc_flyby_SMBH.append(data_file.iloc[0][14])

                        self.freq_flyby_nn.append(data_file.iloc[0][15])
                        self.strain_flyby_nn.append(data_file.iloc[0][16])
                        self.time_flyby_nn.append(data_file.iloc[0][17])
                        self.semi_flyby_nn.append(data_file.iloc[0][18])
                        self.ecc_flyby_nn.append(data_file.iloc[0][19])
                        self.fb_nn_SMBH.append(data_file.iloc[0][20])

                        self.freq_flyby_t.append(data_file.iloc[0][21])
                        self.strain_flyby_t.append(data_file.iloc[0][22])
                        self.time_flyby_t.append(data_file.iloc[0][23])
                        self.semi_flyby_t.append(data_file.iloc[0][24])
                        self.ecc_flyby_t.append(data_file.iloc[0][25])
                        self.fb_t_SMBH.append(data_file.iloc[0][26])

        self.sim_time = np.asarray(self.sim_time, dtype = 'object')
        self.pop = np.asarray(self.pop, dtype = 'object')
        self.mass_SMBH = np.asarray(self.mass_SMBH, dtype = 'object')
        self.mass_parti = np.asarray(self.mass_parti, dtype = 'object')
        self.mass_IMBH = np.asarray(self.mass_IMBH, dtype = 'object')
        self.fb_nn_events = np.asarray(self.fb_nn_events, dtype = 'object')
        self.fb_t_events = np.asarray(self.fb_t_events, dtype = 'object')
        self.tot_events = np.asarray(self.tot_events, dtype = 'object')
        self.close_enc = np.asarray(self.close_enc, dtype = 'object')
        self.freq_flyby_SMBH = np.asarray(self.freq_flyby_SMBH, dtype = 'object')
        self.strain_flyby_SMBH = np.asarray(self.strain_flyby_SMBH, dtype = 'object')
        self.time_flyby_SMBH = np.asarray(self.time_flyby_SMBH, dtype = 'object')
        self.semi_flyby_SMBH = np.asarray(self.semi_flyby_SMBH, dtype = 'object')
        self.ecc_flyby_SMBH = np.asarray(self.ecc_flyby_SMBH, dtype = 'object')
        self.freq_flyby_nn = np.asarray(self.freq_flyby_nn, dtype = 'object')
        self.strain_flyby_nn = np.asarray(self.strain_flyby_nn, dtype = 'object')
        self.time_flyby_nn = np.asarray(self.time_flyby_nn, dtype = 'object')
        self.semi_flyby_nn = np.asarray(self.semi_flyby_nn, dtype = 'object')
        self.ecc_flyby_nn = np.asarray(self.ecc_flyby_nn, dtype = 'object')
        self.fb_nn_SMBH = np.asarray(self.fb_nn_SMBH, dtype = 'object')
        self.freq_flyby_t = np.asarray(self.freq_flyby_t, dtype = 'object')
        self.strain_flyby_t = np.asarray(self.strain_flyby_t, dtype = 'object')
        self.time_flyby_t = np.asarray(self.time_flyby_t, dtype = 'object')
        self.semi_flyby_t = np.asarray(self.semi_flyby_t, dtype = 'object')
        self.ecc_flyby_t = np.asarray(self.ecc_flyby_t, dtype = 'object')
        self.fb_t_SMBH = np.asarray(self.fb_t_SMBH, dtype = 'object')

        with open('figures/gravitational_waves/output/event_rate_'+str(self.integrator[int_idx])+'.txt', 'w') as file:
            tot_events = 0
            tot_IMBH_nn = 0
            tot_IMBH_t = 0
            tot_SMBH_nn = 0
            tot_SMBH_t = 0
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
            
            no_data = np.shape(self.freq_flyby_nn)[0]
            for parti_ in range(no_data):
                sim_time = self.sim_time[parti_]
                parti_nn_filter = np.asarray(self.fb_nn_SMBH[parti_])
                parti_t_filter = np.asarray(self.fb_t_SMBH[parti_])

                no_IMBH_nn = len(np.asarray(self.freq_flyby_nn[parti_])[parti_nn_filter < 0]) / sim_time
                no_SMBH_nn = len(np.asarray(self.freq_flyby_nn[parti_])[parti_nn_filter > 0]) / sim_time
                no_IMBH_t = len(np.asarray(self.freq_flyby_t[parti_])[parti_t_filter < 0]) / sim_time
                no_SMBH_t = len(np.asarray(self.freq_flyby_t[parti_])[parti_t_filter > 0]) / sim_time

                tot_IMBH_nn += no_IMBH_nn / 2
                tot_IMBH_t += no_IMBH_t / 2
                tot_SMBH_nn += no_SMBH_nn
                tot_SMBH_t += no_SMBH_t

                tot_events += len(np.asarray(self.freq_flyby_nn[parti_])[parti_nn_filter < 0]) / sim_time
                tot_events += len(np.asarray(self.freq_flyby_t[parti_])[parti_t_filter < 0]) / sim_time
                tot_events += len(self.freq_flyby_SMBH[parti_]) / sim_time

                freq_temp = np.asarray(self.freq_flyby_nn[parti_])
                strain_temp = np.asarray(self.strain_flyby_nn[parti_])
                fbb_IMBH = np.asarray(self.fb_nn_SMBH[parti_])
                obs_fbnn = freq_temp[(freq_temp > 10**-4) & (fbb_IMBH < 0) & (strain_temp > 10**-27)]
                obs_event_nn += (len(obs_fbnn)/sim_time)
                perc_SMBH_nn += len(obs_fbnn)/len(freq_temp)
                if len(obs_fbnn) > 0:
                    samples_nn += 1
                    no_strain_nn += len(obs_fbnn) / len(freq_temp[(freq_temp > 10**-4) & (fbb_IMBH < 0)])
                
                freq_temp = np.asarray(self.freq_flyby_t[parti_])
                strain_temp = np.asarray(self.strain_flyby_t[parti_])
                fbt_IMBH = np.asarray(self.fb_t_SMBH[parti_])
                obs_fbt = freq_temp[(freq_temp > 10**-4) & (fbt_IMBH < 0) & (strain_temp > 10**-27)]
                obs_event_t += (len(obs_fbt)/sim_time)
                perc_SMBH_t += len(obs_fbt)/len(freq_temp)
                if len(obs_fbt) > 0:
                    samples_t += 1
                    no_strain_t += len(obs_fbt) / len(freq_temp[(freq_temp > 10**-4) & (fbt_IMBH < 0)])

                freq_temp = np.asarray(self.freq_flyby_SMBH[parti_])
                strain_temp = np.asarray(self.strain_flyby_SMBH[parti_])
                obs_fbS = freq_temp[(freq_temp > 10**-4) & (strain_temp > 10**-27)]
                obs_event_SMBH += (len(obs_fbS))/sim_time
                if len(obs_fbS) > 0:
                    samples_SMBH += 1
                    no_strain_S += len(obs_fbS) / len(freq_temp[(freq_temp > 10**-4)])

            obs_event_nn /= max(1, samples_nn)
            obs_event_t /= max(1, samples_t)
            obs_event_SMBH /= max(1, samples_SMBH)
            no_strain_nn /= max(1, no_data)
            no_strain_t /= max(1, no_data)
            no_strain_S /= max(1, no_data)

            fb_nn_events = np.asarray(self.fb_nn_events)
            fb_t_events = np.asarray(self.fb_t_events)
            sim_time = np.asarray(self.sim_time)

            file.write('Stats are for N = 10 and only looks at events with strain > 1e-30 and f > 1e-12 Hz')
            file.write('\nDetectable means f > 1e-4 Hz and strain > 1e-27 when observing a source 1Mpc away\n')
            file.write('\nData for '+str(self.integrator[int_idx]))
            file.write('\nAverage total events per Myr:                       ' + str(tot_events/no_data * 10**6))
            file.write('\nAverage total flyby events per Myr:                 ' + str(np.mean(fb_nn_events/sim_time) * 10**6))
            file.write('\nAverage total tertiary events per Myr:              ' + str(np.mean(fb_t_events/sim_time) * 10**6))
            file.write('\nAverage IMBH-IMBH nn events per Myr:                ' + str(tot_IMBH_nn/max(1, no_data) * 10**6))
            file.write('\nAverage SMBH-IMBH nn events per Myr:                ' + str(tot_SMBH_nn/max(1, no_data) * 10**6))
            file.write('\nAverage IMBH-IMBH tertiary events per Myr:          ' + str(tot_IMBH_t/max(1, no_data) * 10**6))
            file.write('\nAverage SMBH-IMBH tertiary events per Myr:          ' + str(tot_SMBH_t/max(1, no_data) * 10**6))
            file.write('\nAverage detectable nn IMBH-IMBH per Myr:            ' + str(obs_event_nn * 10**6))
            file.write('\nAverage detectable tertiary IMBH-IMBH per Myr:      ' + str(obs_event_t * 10**6))
            file.write('\nAverage detectable IMBH-IMBH per Myr:               ' + str((obs_event_nn + obs_event_t) * 10**6))
            file.write('\nAverage detectable SMBH-IMBH event rate per Myr:    ' + str(obs_event_SMBH * 10**6))
            file.write('\nObservable / Events in freq. range IMBH-IMBH:       ' + str((no_strain_nn+no_strain_t)/2))
            file.write('\nObservable / Events in freq. range SMBH-IMBH:       ' + str(no_strain_S))
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
        self.text_angle = np.degrees(np.arctan((ecc_range[30]-ecc_range[20])/(self.LIGO_semimaj[30]-self.LIGO_semimaj[20])))

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

    def gw_dfreq(self, ecc, m1, m2, semi, chirp_mass, ecc_func):
        """
        Function to take into account the limited LISA observation time, T ~ 5yrs
        Based on equation (6) of Kremer et al. 2019.
        """

        redshift = 0.0228
        nharm = self.gw_harmonic_mode(ecc)
        forb = np.sqrt(constants.G*(m1+m2))/(2*np.pi) * abs(semi)**-1.5 * (redshift+1)**-1

        dfreq = (96*nharm)/(10*np.pi)*(constants.G*chirp_mass)**(5/3)/(constants.c**5) * (2*np.pi * forb)**(11/3) * abs(ecc_func)
        return dfreq

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
        ecc = abs(ecc)

        chirp_mass = (m1*m2)**0.6/(m1+m2)**0.2 * (1+redshift)**-1
        cfactor = 2/(3*np.pi**(4/3)) * (constants.G**(5/3))/(constants.c**3) * (dist*(1+redshift))**-2
        ecc_func = (1+(73/24)*ecc**2+(37/96)*ecc**4)*(1-ecc**2)**-3.5

        nharm = self.gw_harmonic_mode(ecc)
        freq = self.gw_freq(semi, ecc, m1, m2) * (1 | units.Hz)
        dfreq = self.gw_dfreq(ecc, m1, m2, semi, chirp_mass, ecc_func)
        factor = dfreq * (5 * (1 | units.yr))/freq

        strain = min(1, factor) * cfactor * chirp_mass**(5/3) * freq**(-1/3) * (2/nharm)**(2/3) * (self.gfunc(ecc)/ecc_func)
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

    def scatter_hist(self, x1, y1, x2, y2, ax, ax_histf, ax_histh, label1, label2):
        """
        Function to plot the frequency/strain histogram along its scatter plot.
        Use of: https://arxiv.org/pdf/2007.04241.pdf and https://arxiv.org/pdf/1408.0740.pdf 
        
        Inputs:
        x:        The strain arrays
        y:        The frequency arrays
        ax:       axis where the scatter plot is located
        ax_histf: axis where the strain histogram is placed
        ax_histh: axis where the frequency histogram is placed
        len_arr:  Length of SMBH data array
        """

        ######## PLOT {ALL | ONLY SMBH | ONLY NN | ONLY TERT}

        # the scatter plot:
        ax.scatter(np.log10(x1), np.log10(y1), color = 'blueviolet', s = 0.75)
        ax.scatter(np.log10(x2), np.log10(y2), color = 'orange', s = 0.75)
        
        ax_histf.tick_params(axis="x", labelbottom=False)
        ax_histh.tick_params(axis="y", labelleft=False)

        kdef_SMBH = sm.nonparametric.KDEUnivariate(np.log10(x1))
        kdef_SMBH.fit()
        kdef_SMBH.density = (kdef_SMBH.density/max(kdef_SMBH.density))

        kdef_IMBH = sm.nonparametric.KDEUnivariate(np.log10(x2))
        kdef_IMBH.fit()
        kdef_IMBH.density = (kdef_IMBH.density / max(kdef_IMBH.density))

        ax_histf.plot(kdef_SMBH.support, kdef_SMBH.density, color = 'blueviolet', label = label1)
        ax_histf.fill_between(kdef_SMBH.support, kdef_SMBH.density, alpha = 0.35, color = 'blueviolet')
        ax_histf.plot(kdef_IMBH.support, kdef_IMBH.density, color = 'orange', label = label2)
        ax_histf.fill_between(kdef_IMBH.support, kdef_IMBH.density, alpha = 0.35, color = 'orange')
        ax_histf.set_ylim(0, 1.05)
        ax_histf.set_ylabel(r'$\rho/\rho_{\rm{max}}$')
        ax_histf.legend()

        kdeh_SMBH = sm.nonparametric.KDEUnivariate(np.log10(y1))
        kdeh_SMBH.fit()
        kdeh_SMBH.density = (kdeh_SMBH.density / max(kdeh_SMBH.density))

        kdeh_IMBH = sm.nonparametric.KDEUnivariate(np.log10(y2))
        kdeh_IMBH.fit()
        kdeh_IMBH.density = (kdeh_IMBH.density / max(kdeh_IMBH.density))

        ax_histh.plot(kdeh_SMBH.density, kdeh_SMBH.support, color = 'blueviolet')
        ax_histh.fill_between(kdeh_SMBH.density, kdeh_SMBH.support, alpha = 0.35, color = 'blueviolet')
        ax_histh.plot(kdeh_IMBH.density, kdeh_IMBH.support, color = 'orange')
        ax_histh.fill_between(kdeh_IMBH.density, kdeh_IMBH.support, alpha = 0.35, color = 'orange')
        ax_histh.set_xlim(0, 1.05) 
        ax_histh.set_xlabel(r'$\rho/\rho_{\rm{max}}$')

        # LISA
        lisa = li.LISA() 
        x_temp = np.linspace(10**-5, 1, 1000)
        Sn = lisa.Sn(x_temp)

        # SKA
        SKA = np.load(os.path.join(os.path.dirname(__file__), 'SGWBProbe/files/hc_SKA.npz'))
        SKA_freq = SKA['x']
        SKA_hc = SKA['y']
        SKA_strain = SKA_hc**2/SKA_freq

        # BBO 
        BBO = np.load(os.path.join(os.path.dirname(__file__), 'SGWBProbe/files/S_h_BBO_STAR.npz'))
        BBO_freq = BBO['x']
        BBO_strain = BBO['y']

        # muAres 
        Ares = np.load(os.path.join(os.path.dirname(__file__), 'SGWBProbe/files/S_h_muAres_nofgs.npz'))
        Ares_freq = Ares['x']
        Ares_strain = Ares['y']

        ax.plot(np.log10(x_temp), np.log10(np.sqrt(x_temp*Sn)), color = 'slateblue')
        ax.plot(np.log10(BBO_freq), np.log10(np.sqrt(BBO_freq*BBO_strain)), linewidth='1.5', color='blue')
        ax.plot(np.log10(Ares_freq), np.log10(np.sqrt(Ares_freq*Ares_strain)), linewidth='1.5', color='red')
        ax.plot(np.log10(SKA_freq), np.log10(np.sqrt(SKA_freq*SKA_strain)), linewidth='1.5', color='orangered')
        ax.text(-9.25, -15.8, 'SKA', fontsize ='small', rotation = 322, color = 'orangered')
        ax.text(-4.28, -18.2, 'LISA', fontsize ='small', rotation = 309, color = 'slateblue')
        ax.text(-6.13, -19, r'$\mu$Ares', fontsize ='small', rotation = 312, color = 'red')
        ax.text(-1.32, -24, 'BBO', fontsize ='small', rotation = 321, color = 'blue')
        
    def orbital_hist_plotter(self):
        """
        Function which plots all transient events into a histogram.
        Separates events depending on IMBH-IMBH or SMBH-IMBH.
        """

        plot_init = plotter_setup()

        xmin = -10
        xmax = 0
        ymin = -4

        fig = plt.figure(figsize=(13, 12))
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        ax3 = fig.add_subplot(223)
        ax4 = fig.add_subplot(224)
        ax_top = [ax1, ax2]
        ax_bot = [ax3, ax4]
        ax1.set_title('Hermite\nIMBH-IMBH')
        ax2.set_title('GRX\nIMBH-IMBH')
        ax3.set_title('SMBH-IMBH')
        ax4.set_title('SMBH-IMBH')
        extent = [1.1*xmin, 1.1*xmax, 1.1*ymin, 0]
        ecc_l = np.linspace(0, 1, 1000)
        const_semi = [np.log10((1-i)**-1/10)-1 for i in ecc_l]
        const_semi2 = [np.log10((1-i)**-2/10)-2 for i in ecc_l]
        ecc_lx = [np.log10(1-i)+1 for i in ecc_l]

        IMBH_tgw = [[ ], [ ]]
        IMBH_sem = [[ ], [ ]]
        IMBH_ecc = [[ ], [ ]]

        SMBH_tgw = [[ ], [ ]]
        SMBH_sem = [[ ], [ ]]
        SMBH_ecc = [[ ], [ ]]
        
        for int_ in range(2):
            self.combine_data('pop_filt', self.integrator[int_])

            for parti_ in range(len(self.semi_flyby_nn)): #Looping through every individual particle
                for event_ in range(len(self.semi_flyby_nn[parti_])): #Looping through every detected event
                    semi_fb_nn = self.semi_flyby_nn[parti_][event_]
                    if semi_fb_nn < 1 | units.parsec:
                        if np.asarray(self.fb_nn_SMBH[parti_][event_]) < 0:
                            ecc_fb_nn = self.ecc_flyby_nn[parti_][event_]
                            gwtime_nn = self.gw_timescale(semi_fb_nn, ecc_fb_nn, self.mass_parti[parti_], 
                                                          self.mass_parti[parti_]).value_in(units.Myr)

                            IMBH_tgw[int_].append(gwtime_nn)
                            IMBH_sem[int_].append(semi_fb_nn.value_in(units.pc))
                            IMBH_ecc[int_].append(np.log10(1-ecc_fb_nn))
                        else:
                            ecc_fb_nn = self.ecc_flyby_nn[parti_][event_]
                            gwtime_nn = self.gw_timescale(semi_fb_nn, ecc_fb_nn, self.mass_SMBH[parti_], 
                                                            self.mass_parti[parti_]).value_in(units.Myr)

                            SMBH_tgw[int_].append(gwtime_nn)
                            SMBH_sem[int_].append(semi_fb_nn.value_in(units.pc))
                            SMBH_ecc[int_].append(np.log10(1-ecc_fb_nn))

            for parti_ in range(len(self.semi_flyby_t)):
                for event_ in range(len(self.semi_flyby_t[parti_])):
                    semi_fb_t = self.semi_flyby_t[parti_][event_]
                    if semi_fb_t < 1 | units.parsec:
                        if np.asarray(self.fb_t_SMBH[parti_][event_]) < 0 or self.mass_IMBH[parti_] < 10**5 | units.MSun:
                            ecc_fb_t = self.ecc_flyby_t[parti_][event_]
                            gwtime_t = self.gw_timescale(semi_fb_t, ecc_fb_t, self.mass_parti[parti_], 
                                                            self.mass_IMBH[parti_]).value_in(units.Myr)

                            IMBH_tgw[int_].append(gwtime_t)
                            IMBH_sem[int_].append(semi_fb_t.value_in(units.pc))
                            IMBH_ecc[int_].append(np.log10(1-ecc_fb_t))
                        else:
                            ecc_fb_t = self.ecc_flyby_t[parti_][event_]
                            gwtime_t = self.gw_timescale(semi_fb_t, ecc_fb_t, self.mass_parti[parti_], 
                                                            self.mass_SMBH[parti_]).value_in(units.Myr)

                            SMBH_tgw[int_].append(gwtime_t)
                            SMBH_sem[int_].append(semi_fb_t.value_in(units.pc))
                            SMBH_ecc[int_].append(np.log10(1-ecc_fb_t))

            for parti_ in range(len(self.semi_flyby_SMBH)):
                for event_ in range(len(self.semi_flyby_SMBH[parti_])):
                    semi_fb_SMBH = self.semi_flyby_SMBH[parti_][event_]
                    if semi_fb_SMBH < 1 | units.parsec:
                        ecc_fb_SMBH = self.ecc_flyby_SMBH[parti_][event_]
                        gwtime_SMBH = self.gw_timescale(semi_fb_SMBH, ecc_fb_SMBH, self.mass_parti[parti_], 
                                                        self.mass_SMBH[parti_]).value_in(units.Myr)

                        SMBH_tgw[int_].append(gwtime_SMBH)
                        SMBH_sem[int_].append(semi_fb_SMBH.value_in(units.pc))
                        SMBH_ecc[int_].append(np.log10(1-ecc_fb_SMBH))

            IMBH_tgw[int_] = np.asarray(IMBH_tgw[int_])
            IMBH_ecc[int_] = np.asarray(IMBH_ecc[int_])
            IMBH_sem[int_] = np.asarray(IMBH_sem[int_]) 
            SMBH_tgw[int_] = np.asarray(SMBH_tgw[int_])
            SMBH_ecc[int_] = np.asarray(SMBH_ecc[int_])
            SMBH_sem[int_] = np.asarray(SMBH_sem[int_])
        
            ############### PLOTTING OF a vs. (1-e) FOR BIN ##############

            red_mass = (self.mass_parti[0][0]**2)*(2*self.mass_parti[0][0])
            x_arr = np.linspace(-10, 0, 2000)
            const_tgw = [np.log10(1-np.sqrt(1-((256*self.tH*(constants.G**3)/(5*constants.c**5)*red_mass*(10**(i) * (1 | units.pc)) **-4)) **(1/3.5))) for i in x_arr]

            red_mass2 = (self.mass_parti[0][0]*self.mass_SMBH[0][0])*(self.mass_parti[0][0]+self.mass_SMBH[0][0])
            const_tgw2 = [np.log10(1-np.sqrt(1-((256*self.tH*(constants.G**3)/(5*constants.c**5)*red_mass2*(10**(i) * (1 | units.pc)) **-4)) **(1/3.5))) for i in x_arr]

            bin2d_sim, xedg, xedg, image = ax_top[int_].hist2d(np.log10(IMBH_sem[int_]), IMBH_ecc[int_], bins=(150, 150), 
                                                            range=([1.1*xmin, 1.1*xmax], [1.1*ymin, 0]))
            bin2d_sim /= np.max(bin2d_sim)
            contours = ax_top[int_].imshow((bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
            ax_top[int_].scatter(np.log10(IMBH_sem[int_])[np.log10(IMBH_sem[int_]) < -3], 
                                 IMBH_ecc[int_][np.log10(IMBH_sem[int_]) < -3], 
                                 marker = 'X', color = 'white', edgecolor = 'black', s = 5)

            bin2d_sim, xedg, xedg, image = ax_bot[int_].hist2d(np.log10(SMBH_sem[int_]), SMBH_ecc[int_], bins=(150, 150), 
                                                            range=([1.1*xmin, 1.1*xmax], [1.1*ymin, 0]))
            bin2d_sim /= np.max(bin2d_sim)
            contours = ax_bot[int_].imshow((bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
            ax_bot[int_].scatter(np.log10(SMBH_sem[int_])[np.log10(SMBH_sem[int_]) < -3], 
                                 SMBH_ecc[int_][np.log10(SMBH_sem[int_]) < -3], 
                                 marker = 'X', color = 'white', edgecolor = 'black')
            
            for ax_ in [ax_bot[int_], ax_top[int_]]:
                ax_.set_ylabel(r'$\log_{10}(1-e)$')
                ax_.set_xlabel(r'$\log_{10} a$ [pc]')
                ax_.set_xlim(xmin, xmax)
                ax_.set_ylim(ymin, 0)
                plot_init.tickers(ax_, 'hist')

            self.forecast_interferometer(ax_top[int_], self.mass_parti[0][0], self.mass_IMBH[0][0])
            ax_top[int_].plot(x_arr, const_tgw, color = 'white')
            ax_top[int_].text(-9.5, -3, r'aLIGO ($f_{\rm{peak}} = 200$ Hz)', verticalalignment = 'center', fontsize ='small', rotation=self.text_angle-10, color = 'white')
            ax_top[int_].text(-6.6, -3, r'LISA ($f_{\rm{peak}} = 10^{-2}$ Hz)', verticalalignment = 'center', fontsize ='small', rotation=self.text_angle-10, color = 'white')
            ax_top[int_].text(-2.55, -3, r'$t_{\rm{GW}} > t_H$', verticalalignment = 'center', fontsize ='small', rotation=self.text_angle+2, color = 'white')
            ax_top[int_].text(-2.93, -3, r'$t_{\rm{GW}} < t_H$', verticalalignment = 'center', fontsize ='small', rotation=self.text_angle+2, color = 'white')

            self.forecast_interferometer(ax_bot[int_], self.mass_parti[0][0], self.mass_SMBH[0][0])
            ax_bot[int_].plot(x_arr, const_tgw2, color = 'white')
            ax_bot[int_].text(-8.4, -3, r'aLIGO ($f_{\rm{peak}} = 200$ Hz)', verticalalignment = 'center', fontsize ='small', rotation=self.text_angle-10, color = 'white')
            ax_bot[int_].text(-5.5, -3, r'LISA ($f_{\rm{peak}} = 10^{-2}$ Hz)', verticalalignment = 'center', fontsize ='small', rotation=self.text_angle-10, color = 'white')
            ax_bot[int_].text(-0.85, -3, r'$t_{\rm{GW}} > t_H$', verticalalignment = 'center', fontsize ='small', rotation=self.text_angle+4, color = 'white')
            ax_bot[int_].text(-1.23, -3, r'$t_{\rm{GW}} < t_H$', verticalalignment = 'center', fontsize ='small', rotation=self.text_angle+4, color = 'white')
        #ax_top[1].plot(ecc_lx, const_semi)
        #ax_top[0].plot(ecc_lx, const_semi)
        #ax_top[0].plot(ecc_lx, const_semi2, color = 'white')
        plt.savefig('figures/gravitational_waves/ecc_semi_bins_IMBH_histogram_t.png', dpi=500, bbox_inches='tight')
        plt.clf()

        with open('figures/gravitational_waves/output/GW_merger_time.txt', 'w') as file:
            for int_ in range(2):
                file.write('\nData for '+str(self.integrator[int_]))
                file.write('\nAverage avg. GW timescales for IMBH-IMBH:   ' + str(np.mean(IMBH_tgw[int_])) + ' Myr')
                file.write('\nFive smallest GW timescales for IMBH-IMBH:  ' + str(np.sort(IMBH_tgw[int_])[:5]) + ' Myr')
                file.write('\nAverage avg. GW timescales for SMBH-IMBH:   ' + str(np.mean(SMBH_tgw[int_])) + ' Myr')
                file.write('\nFive smallest GW timescales for SMBH-IMBH:  ' + str(np.sort(SMBH_tgw[int_])[:5]) + ' Myr')
                file.write('\n========================================================================')

    def Nenc_tgw_plotter(self):
        """
        Function which manipulates data to extract:
        - Minimum and average tgw w.r.t SMBH for each particle
        - Minimum and average semi-major axis w.r.t SMBH for each particle
        - Minimum and average eccentricity w.r.t SMBH for each particle
        It plots a histogram to convey the evolution of this binary over the course of the simulation.
        It also plots a vs. (1-e) in both a histogram and scatter of the minimum to compare with LISA and LIGO.
        """

        plot_init = plotter_setup()
        
        fig = plt.figure(figsize=(12, 6))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax1.set_title('Hermite')
        ax2.set_title('GRX')

        tgw_SMBH = [[ ], [ ]]
        Nenc = [[ ], [ ]]
        for int_ in range(2):
            self.combine_data('pop_filt', self.integrator[int_])
            for i in range(len(self.semi_flyby_SMBH)):
                tgw_SMBH_avg = []
                for parti_ in range(len(self.semi_flyby_SMBH[i])):
                    grav_avg_time = self.gw_timescale(self.semi_flyby_SMBH[i][parti_], self.ecc_flyby_SMBH[i][parti_], 
                                                    self.mass_SMBH[i], self.mass_SMBH[i]).value_in(units.Myr)
                    tgw_SMBH_avg.append(grav_avg_time)
                data_pts = len(tgw_SMBH_avg)
                tgw_SMBH[int_].append(np.nanmean(np.asarray(tgw_SMBH_avg[round(0.9*data_pts):])))
            tgw_SMBH[int_] = np.asarray(tgw_SMBH[int_])

            for parti_ in range(len(self.close_enc)): #Looping through every individual particle
                Nenc[int_].append(self.close_enc[parti_])

            xmax = (max(Nenc[0]))
            ymin = (np.nanmin(np.log10(tgw_SMBH[0])))
            ymax = (np.nanmax(np.log10(tgw_SMBH[0])))
            iter = -1
            for ax_ in [ax1, ax2]:
                iter += 1
                for int_ in range(2):
                    bin2d_sim, xedg, xedg, image = ax_.hist2d(Nenc[int_], np.log10(tgw_SMBH[int_]), 
                                                        bins=(100,100), range=([0, 1.1*xmax], [0.9*ymin, 1.1*ymax]))
                bin2d_sim /= np.max(bin2d_sim)
                extent = [0, 1.1*xmax, 0.9*ymin, 1.1*ymax]
                contours = ax_.imshow((bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
                ax_.set_xlabel(r'$N_{\rm{enc}}$')
                ax_.set_ylabel(r'$\log_{10}t_{\rm{GW, min}}$ [Myr]')
                ax_.set_xlim(0, 1.1*xmax)
                plot_init.tickers(ax_, 'hist')
        plt.savefig('figures/gravitational_waves/Nenc_tgw_SMBH_evol_histogram.pdf', dpi=300, bbox_inches='tight')
        plt.clf()

    def strain_freq_plotter(self):
        """
        Function which plots the amplitude histogram for all sets of particle.
        """
        
        plot_ini = plotter_setup()
        self.combine_data('integrator', self.integrator[1])

        IMBH_strain = [ ]
        IMBH_strain_nn = [ ]
        IMBH_strain_t = [ ]
        IMBH_freq = [ ]
        IMBH_freq_nn = [ ]
        IMBH_freq_t = [ ]

        SMBH_strain = [ ]
        SMBH_strain_nn = [ ]
        SMBH_strain_t = [ ]
        SMBH_freq = [ ]
        SMBH_freq_nn = [ ]
        SMBH_freq_t = [ ]
        
        for parti_ in range(len(self.semi_flyby_nn)): #Looping through every individual particle
            for event_ in range(len(self.semi_flyby_nn[parti_])): #Looping through every detected event
                semi_fb_nn = self.semi_flyby_nn[parti_][event_]
                if semi_fb_nn < 1 | units.parsec and self.strain_flyby_nn[parti_][event_] > 0:
                    if np.asarray(self.fb_nn_SMBH[parti_][event_]) < 0:
                        IMBH_strain.append(self.strain_flyby_nn[parti_][event_])
                        IMBH_strain_nn.append(self.strain_flyby_nn[parti_][event_])
                        IMBH_freq.append(self.freq_flyby_nn[parti_][event_])
                        IMBH_freq_nn.append(self.freq_flyby_nn[parti_][event_])
                    else:
                        SMBH_strain.append(self.strain_flyby_nn[parti_][event_])
                        SMBH_strain_nn.append(self.strain_flyby_nn[parti_][event_])
                        SMBH_freq.append(self.freq_flyby_nn[parti_][event_])
                        SMBH_freq_nn.append(self.freq_flyby_nn[parti_][event_])

        for parti_ in range(len(self.semi_flyby_t)):
            for event_ in range(len(self.semi_flyby_t[parti_])):
                semi_fb_t = self.semi_flyby_t[parti_][event_]
                if semi_fb_t < 1 | units.parsec and self.strain_flyby_t[parti_][event_] > 0:
                    if np.asarray(self.fb_t_SMBH[parti_][event_]) < 0:
                        IMBH_strain.append(self.strain_flyby_t[parti_][event_])
                        IMBH_strain_t.append(self.strain_flyby_t[parti_][event_])
                        IMBH_freq.append(self.freq_flyby_t[parti_][event_])
                        IMBH_freq_t.append(self.freq_flyby_t[parti_][event_])
                    else:
                        SMBH_strain.append(self.strain_flyby_t[parti_][event_])
                        SMBH_strain_t.append(self.strain_flyby_t[parti_][event_])
                        SMBH_freq.append(self.freq_flyby_t[parti_][event_])
                        SMBH_freq_t.append(self.freq_flyby_t[parti_][event_])

        for parti_ in range(len(self.semi_flyby_SMBH)):
            for event_ in range(len(self.semi_flyby_SMBH[parti_])):
                semi_fb_SMBH = self.semi_flyby_SMBH[parti_][event_]
                if semi_fb_SMBH < 1 | units.parsec:
                    SMBH_strain.append(self.strain_flyby_SMBH[parti_][event_])
                    SMBH_freq.append(self.freq_flyby_SMBH[parti_][event_])

        fig = plt.figure(figsize=(8, 6))
        gs = fig.add_gridspec(2, 2,  width_ratios=(4, 2), height_ratios=(2, 4),
                              left=0.1, right=0.9, bottom=0.1, top=0.9,
                              wspace=0.05, hspace=0.05)
        ax = fig.add_subplot(gs[1, 0])
        ax1 = fig.add_subplot(gs[0, 0], sharex=ax)
        ax2 = fig.add_subplot(gs[1, 1], sharey=ax)
        
        self.scatter_hist(SMBH_freq, SMBH_strain,
                          IMBH_freq, IMBH_strain, 
                          ax, ax1, ax2, 'SMBH-IMBH', 'IMBH-IMBH')
        ax.set_xlabel(r'$\log_{10}f$ [Hz]')
        ax.set_ylabel(r'$\log_{10}h$')
        ax1.set_title(str(self.integrator[1]))
        plot_ini.tickers(ax, 'plot')
        plot_ini.tickers(ax1, 'plot')
        plot_ini.tickers(ax2, 'plot')
        ax.set_ylim(-30, -12.2)
        ax.set_xlim(-12.5, 0.1)
        plt.savefig('figures/gravitational_waves/'+str(self.integrator[1])+'GW_freq_strain_maximise_diagram.png', dpi = 500, bbox_inches='tight')
        plt.clf()

        SMBH_mass = self.mass_SMBH[0][0]
        IMBH_mass = self.mass_IMBH[0][0]
        sem_SMBH = self.coll_radius(SMBH_mass)
        sem_IMBH = self.coll_radius(IMBH_mass)

        dist_z = 7015 | units.Mpc
        max_freq_SMBH = (np.pi)**-1 * np.sqrt(constants.G*(SMBH_mass + IMBH_mass)/sem_SMBH**3)  # Eqn. (3) Matsubayashi 2003
        max_strain_SMBH = np.sqrt(32/5) * (np.pi**2*constants.G**5)**(1/3) * constants.c**(-4)*(SMBH_mass * IMBH_mass)*(SMBH_mass + IMBH_mass)**(-1/3)*max_freq_SMBH**(2/3)*dist_z**-1
        max_strain_SMBH = max_strain_SMBH.value_in(units.kg**1.6653345369377348e-16 * units.s**3.3306690738754696e-16)
        max_freq_IMBH = (np.pi)**-1 * np.sqrt(constants.G*(IMBH_mass + IMBH_mass)/sem_IMBH**3)
        max_strain_IMBH = np.sqrt(32/5) * (np.pi**2*constants.G**5)**(1/3) * constants.c**(-4)*(IMBH_mass * IMBH_mass)*(IMBH_mass + IMBH_mass)**(-1/3)*max_freq_IMBH**(2/3)*dist_z**-1
        max_strain_IMBH = max_strain_IMBH.value_in(units.kg**1.6653345369377348e-16 * units.s**3.3306690738754696e-16)

        with open('figures/gravitational_waves/output/Binaries_redshift_freq_strain.txt', 'w') as file:
            IMBH_strain = np.asarray(IMBH_strain)
            SMBH_strain = np.asarray(SMBH_strain)
            SMBH_freq = np.asarray(SMBH_freq)
            file.write('\nData for '+str(self.integrator[1]))
            file.write('\nFor eccentric events @ z ~ 3.5, IMBH-IMBH min. strain is: ' + str(np.nanmin(IMBH_strain[IMBH_strain > 0]) * (1/7015)))
            file.write('\nFor eccentric events @ z ~ 3.5, IMBH-IMBH max. strain is: ' + str(np.nanmax(IMBH_strain[IMBH_strain > 0]) * (1/7015)))
            file.write('\nFor eccentric events @ z ~ 3.5, SMBH-IMBH min. strain is: ' + str(np.nanmin(SMBH_strain[SMBH_strain > 0]) * (1/7015)))
            file.write('\nFor eccentric events @ z ~ 3.5, SMBH-IMBH max. strain is: ' + str(np.nanmax(SMBH_strain[SMBH_strain > 0]) * (1/7015)))
            file.write('\nFor eccentric events @ z ~ 3.5, IMBH-IMBH min. freq is:   ' + str(np.nanmin(IMBH_freq) * (1/7015)))
            file.write('\nFor eccentric events @ z ~ 3.5, IMBH-IMBH max. freq is:   ' + str(np.nanmax(IMBH_freq) * (1/7015)))
            file.write('\nFor eccentric events @ z ~ 3.5, SMBH-IMBH min. freq is:   ' + str(np.nanmin(SMBH_freq[SMBH_freq > 0]) * (1/7015)))
            file.write('\nFor eccentric events @ z ~ 3.5, SMBH-IMBH max. freq is:   ' + str(np.nanmax(SMBH_freq[SMBH_freq > 0]) * (1/7015)))
            file.write('\nInspiral IMBH-IMBH freq is: ' + str(max_freq_IMBH.value_in(units.Hz)) + ' with strain at 1 Mpc: ' + str(max_strain_IMBH))
            file.write('\nInspiral SMBH-IMBH freq is: ' + str(max_freq_SMBH.value_in(units.Hz)) + ' with strain at 1 Mpc: ' + str(max_strain_SMBH))
            file.write('\n========================================================================')

    def transient_events(self):

        plot_init = plotter_setup() 
        self.combine_data('integrator', self.integrator[1])

        time_IMBH_arr = [ ]
        time_SMBH_arr = [ ]
        freq_IMBH_arr = [ ]
        freq_SMBH_arr = [ ]

        for parti_ in range(len(self.time_flyby_nn)):
            if self.pop[parti_] == 10:
                for event_ in range(len(self.time_flyby_nn[parti_])):
                    if self.freq_flyby_nn[parti_][event_] > 10**-5:
                        if self.fb_nn_SMBH[parti_][event_] < 0:
                            time_IMBH_arr.append(self.time_flyby_nn[parti_][event_])
                            freq_IMBH_arr.append(self.freq_flyby_nn[parti_][event_])
                        else:
                            time_SMBH_arr.append(self.time_flyby_nn[parti_][event_])
                            freq_SMBH_arr.append(self.freq_flyby_nn[parti_][event_])

        for parti_ in range(len(self.time_flyby_t)):
            if self.pop[parti_] == 10:
                for event_ in range(len(self.time_flyby_t[parti_])):
                    if self.freq_flyby_t[parti_][event_] > 10**-5:
                        if self.fb_t_SMBH[parti_][event_] < 0:
                            time_IMBH_arr.append(self.time_flyby_t[parti_][event_])
                            freq_IMBH_arr.append(self.freq_flyby_t[parti_][event_])
                        else:
                            time_SMBH_arr.append(self.time_flyby_t[parti_][event_])
                            freq_SMBH_arr.append(self.freq_flyby_t[parti_][event_])

        for parti_ in range(len(self.time_flyby_t)):
            if self.pop[parti_] == 10:
                for event_ in range(len(self.time_flyby_t[parti_])):
                    if self.freq_flyby_SMBH[parti_][event_] > 10**-5:
                        for event_ in range(len(self.time_flyby_SMBH[parti_])):
                            time_SMBH_arr.append(self.time_flyby_SMBH[parti_][event_])
                            freq_SMBH_arr.append(self.freq_flyby_SMBH[parti_][event_])

        fig = plt.figure(figsize=(16, 14))
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        plot_title = ['SMBH-IMBH', 'IMBH-IMBH']
        colours = ['red', 'blue']
        ax = [ax1, ax2]
        iter = 0
        for ax_ in ax:
            plot_init.tickers(ax_, 'plot')
            ax_.set_xlabel('Time [Myr]')
            ax_.set_ylabel(r' $\log_{10} f_{\rm{GW}}$ [Hz]')
            ax_.set_title(plot_title[iter])
            ax_.set_ylim(-5, 0)
            iter += 1

        ax1.scatter(time_SMBH_arr,  np.log10(freq_SMBH_arr), 
                    edgecolors = 'black', color = colours[1])
        ax2.scatter(time_IMBH_arr, np.log10(freq_IMBH_arr), 
                    edgecolors = 'black', color = colours[1])

        plt.savefig('figures/gravitational_waves/events_time.png', dpi = 300, bbox_inches='tight')

    def strain_plotters(self):

        plot_init = plotter_setup() 

        SMBH_mass = self.mass_SMBH[0][0][0]
        IMBH_mass = self.mass_SMBH[0][1][0]
        ecc = [0.1, 0.5, 0.9]
        semi_major = np.linspace(0.001, 0.4, 5000)
        colors = ['black', 'red', 'blue']

        fig, ax = plt.subplots()
        plot_init.tickers(ax, 'plot')
        ax.set_xlabel(r'$\log_{10} a$ [pc]')
        ax.set_ylabel(r'$\log_{10} h$')
        for ecc_ in range(len(ecc)):
            strains = [self.gw_strain(i * (1 | units.pc), ecc[ecc_], IMBH_mass, SMBH_mass) for i in semi_major]
            ax.plot(np.log10(semi_major), np.log10(strains), label = r'$e = $'+str(ecc[ecc_]), color = colors[ecc_])
        ax.legend()
        plt.savefig('figures/gravitational_waves/strain_dependence.pdf', dpi = 300, bbox_inches='tight')


"""#print('...tGW_plotters...')
cst = gw_calcs()
#cst.new_data_extractor()
#cst.orbital_hist_plotter()
#cst.Nenc_tgw_plotter()
#cst.strain_freq_plotter()
cst.transient_events()"""