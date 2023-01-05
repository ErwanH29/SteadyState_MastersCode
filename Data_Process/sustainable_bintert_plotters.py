from amuse.lab import *
from file_logistics import *
from tGW_plotters import *
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import warnings
import pickle as pkl

class sustainable_sys(object):
    """
    Class which extracts ALL information of particle_trajectory files in an memory-optimised manner
    and appends the required data values into arrays.
    """

    def __init__(self):
        warnings.filterwarnings("ignore", category=RuntimeWarning) 

    def new_data_extractor(self):
        """
        Script to extract data from recently simulated runs
        """
        
        GW_calcs = gw_calcs()

        print('!!!!!! WARNING THIS WILL TAKE A WHILE !!!!!!!')
        filenameH = glob.glob(os.path.join('/media/erwanh/Elements/Hermite/particle_trajectory_temp/*'))
        filenameGRX = glob.glob('/media/erwanh/Elements/GRX/particle_trajectory_temp/*')
        filename = [natsort.natsorted(filenameH), natsort.natsorted(filenameGRX)]
        ints = ['Hermite', 'GRX']
        count = 0
        for int_ in range(2):
            for file_ in range(len(filename[int_])):
                with open(filename[int_][file_], 'rb') as input_file:
                    print('Reading file', file_, ':', input_file)
                    file_size = os.path.getsize(filename[int_][file_])
                    if file_size < 2.8e9:
                        count += 1
                        data = pkl.load(input_file)
                        pop = 5*round(0.2*np.shape(data)[0])

                        if pop > 5:
                            for parti_ in range(np.shape(data)[0]):
                                integrator = [ ]
                                if int_ == 0:
                                    integrator.append('Hermite')
                                else:
                                    integrator.append('GRX')

                                pop_bin = [ ]
                                pop_ter = [ ]

                                semi_NN_avg = [ ]
                                semi_NN_min = [ ]
                                semi_t_avg = [ ]
                                semi_t_min = [ ]

                                sys_bin = [ ]
                                bsys_time = [ ]
                                bform_time = [ ]
                                bin_key = [ ]
                                sys_ter = [ ]
                                tsys_time = [ ]
                                tform_time = [ ]
                                ter_key = [ ]

                                GW_freqbin = [ ]
                                GW_strainbin = [ ]
                                GW_timeb = [ ]
                                GW_bmass = [ ]
                                GW_freqter = [ ]
                                GW_strainter = [ ]
                                GW_timet = [ ]
                                GW_tmass = [ ]

                                hard_bin = [ ]
                                hard_ter = [ ]
                                dedt = [ ]
                                dadt = [ ]

                                bin_sys = False
                                ter_sys = False
                                semi_nn_avg_temp = []
                                semi_t_avg_temp = []
                                temp_dedt = [] 
                                temp_dadt = []
                                bform_time.append(-5)
                                tform_time.append(-5)

                                if parti_ != 0:
                                    temp_dedt.append((data.iloc[parti_][-2][8][0] - data.iloc[parti_][2][8][0])/(np.shape(data)[1]-3))
                                    temp_dadt.append((data.iloc[parti_][-2][7][0]**-1 - data.iloc[parti_][2][7][0]**-1).value_in((units.pc)**-1)/(np.shape(data)[1]-3))
                                    for col_ in range(np.shape(data)[1]-1):
                                        nn_semi = abs(data.iloc[parti_][col_][7][1])
                                        nn_ecc = data.iloc[parti_][col_][8][1]
                                        mass1 = data.iloc[parti_][0][1]

                                        for part_ in range(np.shape(data)[0]):
                                            if data.iloc[part_][0][0] == data.iloc[parti_][col_][6][1]:
                                                mass2 = data.iloc[part_][0][1]
                                                mass = max(mass1, mass2)
                                                hard = False
                                                bin = False
                                                if nn_semi < (constants.G*mass)/(4*(150000*(1 | units.ms))**2) and nn_ecc < 1:   #Hard binary conditions (Quinlan 1996b)
                                                    hard = True
                                                    bin = True
                                                    hard_bin.append(1)
                                                    
                                                if not (hard) and nn_semi < (constants.G*mass)/(4*(15000*(1 | units.ms))**2) and abs(nn_ecc) < 1:  #Value chosen for a < 1000AU
                                                    hard_bin.append(-5)
                                                    bin = True

                                                if (bin):
                                                    if not (bin_sys):  #First formation time
                                                        formation_time = col_*1000
                                                        bform_time[-1] = formation_time
                                                        bin_sys = True
                                                    bin_val += 1
                                                    bin_key.append(data.iloc[parti_][col_][6][1])
                                                    bsys_time.append(col_)
                                                    semi_nn_avg_temp.append(nn_semi.value_in(units.pc))
                                                    pop_bin.append(pop)

                                                    strain = GW_calcs.gw_strain(nn_semi, nn_ecc, mass1, mass2)
                                                    freq = GW_calcs.gw_freq(nn_semi, nn_ecc, mass1, mass2)
                                                    GW_time = GW_calcs.gw_timescale(nn_semi, nn_ecc, mass1, mass2)
                                                    GW_strainbin.append(strain)
                                                    GW_freqbin.append(freq)
                                                    GW_timeb.append(float(GW_time.value_in(units.Myr)))
                                                    GW_bmass.append(mass2.value_in(units.MSun))

                                                    semi_outer = abs(data.iloc[parti_][col_][7][2])
                                                    ecc_outer = data.iloc[parti_][col_][8][2]

                                                    #Calculate tertiary. The stability equality is based on Mardling and Aarseth 2001
                                                    if ecc_outer < 1 and semi_outer < 1e-2 | units.parsec:
                                                        for part_ in range(np.shape(data)[0]):
                                                            if data.iloc[part_][0][0] == data.iloc[parti_][col_][6][2]:
                                                                mass_outer = data.iloc[part_][0][1]

                                                        semi_ratio = semi_outer / nn_semi
                                                        equality = 2.8 * ((1+mass_outer/(mass1+mass2))*(1+ecc_outer)/(1-ecc_outer)**0.5)**0.4
                                                        if semi_ratio > equality:
                                                            if not (ter_sys):
                                                                formation_time = col_*1000
                                                                tform_time[-1] = formation_time
                                                                ter_sys = True
                                                            ter_val += 1
                                                            ter_key.append([i for i in data.iloc[parti_][col_][6][2]][0])
                                                            tsys_time.append(col_)
                                                            semi_t_avg_temp.append(semi_outer.value_in(units.pc))
                                                            pop_ter.append(pop)

                                                            GW_time = GW_calcs.gw_timescale(semi_outer, ecc_outer, mass1, mass_outer)
                                                            strain = GW_calcs.gw_strain(semi_outer, ecc_outer, mass1, mass_outer)
                                                            freq = GW_calcs.gw_freq(semi_outer, ecc_outer, mass1, mass_outer)
                                                            GW_strainter.append(strain)
                                                            GW_freqter.append(freq)
                                                            GW_timet.append(float(GW_time.value_in(units.Myr)))
                                                            GW_tmass.append([mass2.value_in(units.MSun), mass_outer.value_in(units.MSun)])

                                                            if (hard):
                                                                hard_ter.append(1)
                                                            else:
                                                                hard_ter.append(-5)

                                        bin_val = 0 
                                        ter_val = 0
                                        bin_sys = 0 
                                        ter_sys = 0

                                        if len(semi_nn_avg_temp) > 0:
                                            semi_NN_avg.append(np.mean(semi_nn_avg_temp))
                                            semi_NN_min.append(np.min(semi_nn_avg_temp))
                                        else:
                                            semi_NN_avg.append(-5)
                                            semi_NN_min.append(-5)
                                        if len(semi_t_avg_temp) > 0:
                                            semi_t_avg.append(np.mean(semi_t_avg_temp))
                                            semi_t_min.append(np.min(semi_t_avg_temp))
                                        else:
                                            semi_t_avg.append(-5)
                                            semi_t_min.append(-5)

                                        if parti_ == (np.shape(data)[0]-1):
                                            if len(bin_key) > 0:
                                                bin_formed = len(np.unique(bin_key)) - 1
                                                bin_sys += bin_formed
                                            if len(ter_key) > 0:
                                                ter_formed = len(np.unique(ter_key)) - 1
                                                ter_sys += ter_formed
                                                
                                    bsys_time.append((len(np.unique(bsys_time)))/(col_))
                                    tsys_time.append((len(np.unique(tsys_time)))/(col_))
                                    sys_bin.append(bin_sys)
                                    sys_ter.append(ter_sys)
                                    dedt.append(np.mean(temp_dedt))
                                    dadt.append(np.mean(temp_dadt))

                                    path = '/home/erwanh/Desktop/SteadyStateBH/Data_Process/data/bin_hier_systems/'
                                    stab_tracker = pd.DataFrame()
                                    df_stabtime = pd.Series({'Integrator': integrator,
                                                            'Population': pop,
                                                            'Binary Pop.': pop_bin,
                                                            '# Binary Sys.': bin_sys,
                                                            'First bin. form': bform_time,
                                                            'Tertiary Pop.': pop_ter,
                                                            '# Tertiary Sys.': ter_sys,
                                                            'First ter. form': tform_time,
                                                            'Hard Bin. Bool': hard_bin,
                                                            'Hard Ter. Bool': hard_ter,
                                                            'Bin. GW freq': GW_freqbin,
                                                            'Bin. GW strain': GW_strainbin,
                                                            'Bin. GW time': GW_timeb,           #Merging time in Myr
                                                            'Bin. GW mass': GW_bmass,           #Constituent mass in MSun
                                                            'Ter. GW freq': GW_freqter,
                                                            'Ter. GW strain': GW_strainter,
                                                            'Ter. GW time': GW_timet,
                                                            'Ter. GW mass': GW_tmass,
                                                            'Average dedt': dedt,
                                                            'Average dadt': dadt,
                                                            'Bin. Semi-major Avg': semi_NN_avg,
                                                            'Bin. Semi-major Min': semi_NN_min,
                                                            'Ter. Semi-major Avg': semi_t_avg,
                                                            'Ter. Semi-major Min': semi_t_min
                                                            })
                                    stab_tracker = stab_tracker.append(df_stabtime, ignore_index = True)
                                    stab_tracker.to_pickle(os.path.join(path, 'IMBH_'+str(ints[int_])+'_system_data_indiv_parti_'+str(count)+'_'+str(parti_)+'_local2.pkl'))

    def array_rewrite(self, arr, new_arr, arr_type, filt):
        """
        Function to rewrite array to manipulatable float format

        Inputs:
        arr:      The original array with the data
        new_arr:  The new, manipulatable, array
        arr_type: String stating whether array is nested or not 
        filt:     Boolean to filter out unwanted values
        """

        if arr_type == 'nested':
            for sublist in arr:
                for item_ in sublist:
                    if (filt):
                        if item_ > 0:
                            new_arr.append(item_)
                    else:
                        new_arr.append(item_)
        else:
            for item_ in arr:
                if (filt):
                    if item_ > 0:
                        new_arr.append(item_)
                else:
                    new_arr.append(item_)
        
        return new_arr
                

    def combine_data(self):
        """
        Function which extracts ALL data and computes
        various quantitative results.
        """

        self.integrator = [[ ], [ ]]
        self.pop = [[ ], [ ]]
        self.pop_bin = [[ ], [ ]]
        self.bin_sys = [[ ], [ ]]
        self.bform_time = [[ ], [ ]]
        self.pop_ter = [[ ], [ ]]
        self.ter_sys = [[ ], [ ]]
        self.tform_time = [[ ], [ ]]
        self.hard_bin = [[ ], [ ]]
        self.hard_ter = [[ ], [ ]]
        self.GW_freqbin = [[ ], [ ]]
        self.GW_strainbin = [[ ], [ ]]
        self.GW_timeb = [[ ], [ ]]
        self.GW_bmass = [[ ], [ ]]
        self.GW_freqter = [[ ], [ ]]
        self.GW_strainter = [[ ], [ ]]
        self.GW_timet = [[ ], [ ]]
        self.GW_tmass = [[ ], [ ]]
        self.dedt = [[ ], [ ]]
        self.dadt = [[ ], [ ]]
        self.semi_NN_avg = [[ ], [ ]]
        self.semi_NN_min = [[ ], [ ]]
        self.semi_t_avg = [[ ], [ ]]
        self.semi_t_min = [[ ], [ ]]

        system_data = natsort.natsorted(glob.glob('data/bin_hier_systems/*'))
        for file_ in range(len(system_data)):
            with open(system_data[file_], 'rb') as input_file:
                data_file = pkl.load(input_file)
                if data_file.iloc[0][0][0] == 'Hermite':
                    int_ = 0
                else:
                    int_ = 1
                self.integrator[int_].append(data_file.iloc[0][0][0])
                self.pop[int_].append(int(data_file.iloc[0][1]))
                self.pop_bin[int_].append(data_file.iloc[0][2])
                self.bin_sys[int_].append(data_file.iloc[0][3])
                self.bform_time[int_].append(data_file.iloc[0][4])
                self.pop_ter[int_].append(data_file.iloc[0][5])
                self.ter_sys[int_].append(data_file.iloc[0][6])
                self.tform_time[int_].append(data_file.iloc[0][7])
                self.hard_bin[int_].append(data_file.iloc[0][8])
                self.hard_ter[int_].append(data_file.iloc[0][9])
                self.GW_freqbin[int_].append(data_file.iloc[0][10])
                self.GW_strainbin[int_].append(data_file.iloc[0][11])
                self.GW_timeb[int_].append(data_file.iloc[0][12])
                self.GW_bmass[int_].append(data_file.iloc[0][13])
                self.GW_freqter[int_].append(data_file.iloc[0][14])
                self.GW_strainter[int_].append(data_file.iloc[0][15])
                self.GW_timet[int_].append(data_file.iloc[0][16])
                self.GW_tmass[int_].append(data_file.iloc[0][17])
                self.dedt[int_].append(data_file.iloc[0][18])
                self.dadt[int_].append(data_file.iloc[0][19])
                self.semi_NN_avg[int_].append(data_file.iloc[0][20])
                self.semi_NN_min[int_].append(data_file.iloc[0][21])
                self.semi_t_avg[int_].append(data_file.iloc[0][22])
                self.semi_t_min[int_].append(data_file.iloc[0][23])
        
        with open('figures/binary_hierarchical/output/system_summary.txt', 'w') as file:
            integrator = ['Hermite', 'GRX']
            for int_ in range(2):
                bform_stat = [ ]
                tform_stat = [ ]
                GWtime_stat = [ ]
                minGWtime = [ ]
                semi_NN_astat = [ ]
                semi_NN_mstat = [ ]
                semi_t_astat = [ ]
                semi_t_mstat = [ ]

                pop_arr = np.unique(self.pop[int_])
                file.write('\n\nData for '+str(integrator[int_]+' in pc'))
                for pop_ in pop_arr:
                    if pop_ < 50:
                        idx = np.argwhere(self.pop[int_] == pop_).flatten()
                        bform_time = [ ]
                        tform_time = [ ]
                        GWb_time = [ ]
                        semi_NN_avg = [ ]
                        semi_NN_min = [ ]
                        semi_t_avg = [ ]
                        semi_t_min = [ ]
                        for data_ in idx:
                            for item_ in self.bform_time[int_][data_]:
                                if item_ > 0:
                                    bform_time.append(item_)
                            for item_ in self.tform_time[int_][data_]:
                                if item_ > 0:
                                    tform_time.append(item_)
                            for item_ in self.GW_timeb[int_][data_]:
                                for k_ in self.hard_bin[int_][data_]:
                                    if k_ > 0 and item_ > 0:
                                        GWb_time.append(float(item_))
                                else:
                                    GWb_time.append(-5)

                            for item_ in self.semi_NN_avg[int_][data_]:
                                if item_ > 0:
                                    semi_NN_avg.append(item_)
                            for item_ in self.semi_NN_min[int_][data_]:
                                if item_ > 0:
                                    semi_NN_min.append(item_)
                            for item_ in self.semi_t_avg[int_][data_]:
                                if item_ > 0:
                                    semi_t_avg.append(item_)
                            for item_ in self.semi_t_min[int_][data_]:
                                if item_ > 0:
                                    semi_t_min.append(item_)
                    GWb_time = [float(i) for i in GWb_time if i > 0]
                    
                    bform_stat.append('{:.7f}'.format(np.mean(bform_time)))
                    tform_stat.append('{:.7f}'.format(np.mean(tform_time)))
                    GWtime_stat.append('{:.7f}'.format(np.mean(GWb_time)))
                    minGWtime.append(np.sort(GWb_time)[:5])
                    semi_NN_astat.append('{:.7f}'.format(np.mean(semi_NN_avg)))
                    semi_NN_mstat.append('{:.7f}'.format(np.mean(semi_NN_min)))
                    semi_t_astat.append('{:.7f}'.format(np.mean(semi_t_avg)))
                    semi_t_mstat.append('{:.7f}'.format(np.mean(semi_t_min)))

                file.write('\nAverage binary formation time [yrs]:             '+str(pop_arr)+' : '+str(bform_stat))
                file.write('\nAverage tertiary formation time [yrs]:           '+str(pop_arr)+' : '+str(tform_stat))
                file.write('\nAverage binary GW timescale [Myr]:               '+str(pop_arr)+' : '+str(GWtime_stat))
                file.write('\nMinimum binary GW timescale [Myr]:               '+str(pop_arr)+' : '+str(minGWtime))
                file.write('\nAverage binary semi-major axis per population:   '+str(pop_arr)+' : '+str(semi_NN_astat))
                file.write('\nMinimum binary semi-major axis per population:   '+str(pop_arr)+' : '+str(semi_NN_mstat))
                file.write('\nAverage tertiary semi-major axis per population: '+str(pop_arr)+' : '+str(semi_t_astat))
                file.write('\nMinimum tertiary semi-major axis per population: '+str(pop_arr)+' : '+str(semi_t_mstat))

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
        
        semi_avNN = [ ]
        semi_avter = [ ]
        bin_pop = [ ]
        GWb_mergertime = [ ]
        GWt_mergertime = [ ]
        GW_bmass = [ ]
        GW_tmass = [ ]
        
        semi_avNN = self.array_rewrite(self.semi_NN_avg[int_], semi_avNN, 'nested', True)
        semi_avter = self.array_rewrite(self.semi_t_avg[int_], semi_avter, 'nested', True)
        bin_pop = self.array_rewrite(self.pop_bin[int_], bin_pop, 'nested', False)
        GWb_mergertime = self.array_rewrite(self.GW_timeb[int_], GWb_mergertime, 'nested', False)
        GWt_mergertime = self.array_rewrite(self.GW_timet[int_], GWt_mergertime, 'nested', False)
        GW_bmass = self.array_rewrite(self.GW_bmass[int_], GW_bmass, 'nested', False)
        GW_tmass = self.array_rewrite(self.GW_tmass[int_], GW_tmass, 'nested', False)

        GWb_mergertime = np.asarray(GWb_mergertime)
        GWt_mergertime = np.asarray(GWt_mergertime)
        GW_bmass = np.asarray(GW_bmass)
        GW_tmass = np.asarray(GW_tmass)

        GW_tmass_in = [ ]
        GW_tmass_out = [ ]
        for i in range(len(GW_tmass)):
            if i%2 == 0:
                GW_tmass_in.append(GW_tmass[i])
            else:
                GW_tmass_out.append(GW_tmass[i])
        GW_tmass_in = np.asarray(GW_tmass_in)
        GW_tmass_out = np.asarray(GW_tmass_out)

        bin_sys = [ ]
        ter_sys = [ ]
        bform_time = [ ]
        tform_time = [ ]

        bin_sys = self.array_rewrite(self.bin_sys[int_], bin_sys, 'not', False)
        ter_sys = self.array_rewrite(self.ter_sys[int_], ter_sys, 'not', False)

        for item_ in self.bform_time[int_][0]:
            if item_ > 0:
                bform_time.append(item_)
        for item_ in self.tform_time[int_][0]:
            if item_ > 0:
                tform_time.append(item_)

        for pop_ in ini_pop:
            iter += 1
            idx2 = np.where(bin_pop == pop_)[0]
            
            bin_formed[iter] = np.mean(bin_sys)
            ter_formed[iter] = np.mean(ter_sys)
            bsys_time[iter] = np.mean(bform_time)
            tsys_time[iter] = np.mean(tform_time)
            bform_ini[iter] = len(np.asarray(bsys_time)[(bsys_time) == 0])
            tform_ini[iter] = len(np.asarray(tsys_time)[(tsys_time) == 0])
        
            if len(semi_avNN) > 0:
                semi_major_bin[iter] = np.mean(semi_avNN)
            if len(semi_avter) > 0:
                semi_major_ter[iter] = np.mean(semi_avter)

        with open('figures/binary_hierarchical/output/'+str(integrator[int_])+'bin_ter_systems.txt', 'w') as file:
            file.write(str(integrator[int_])+' first binary avg. formation time')
            for pop_ in range(len(bform_time)):
                file.write('\nPopulation: '+str(ini_pop[pop_])+':         '+str(bform_time[pop_]/10**3)+' kyr')
            file.write('\n# Binary systems at initialisation:             '+str(bform_ini[pop_]))
            file.write('\n5 shortest binary GW merger time_scales [Myr]:  ')

            for i in range(5):
                file.write('{:.2f}'.format(np.sort(GWb_mergertime[GWb_mergertime > 0])[i])+',   ')
            file.write('\nFraction of mergers within Hubble time:         '+str(np.shape(GWb_mergertime[GWb_mergertime < tH.value_in(units.Myr)])[0]) +'/'+str(np.shape(GWb_mergertime)[0]))
            file.write('\nFraction of IMBH-IMBH binaries:                 '+str(np.shape(GW_bmass[GW_bmass < 10**6])[0]) +'/'+str(np.shape(GW_bmass)[0]))

            file.write('\n\nfirst tertiary avg. formation time')
            for pop_ in range(len(tform_time)):
                file.write('\nPopulation: '+str(ini_pop[pop_])+':          '+str(tform_time[pop_]/10**3)+' kyr')
            file.write('\n# Tertiary systems at initialisation:            '+str(tform_ini[pop_]))
            file.write('\n5 shortest tertiary GW merger time scales [Myr]: ')
            if len(GWt_mergertime) > 5:
                for i in range(5):
                    file.write('{:.2f}'.format(np.sort(GWt_mergertime[GWt_mergertime > 0])[i])+',   ')

            file.write('\nFraction of mergers within Hubble time:          '+str(np.shape(GWt_mergertime[GWt_mergertime < tH.value_in(units.Myr)])[0]) +'/'+str(np.shape(GWt_mergertime)[0]))
            file.write('\nFraction of (IMBH-IMBH)-SMBH tertiaries:         '+str(np.shape(GW_tmass[(GW_tmass_in < 10**6) & (GW_tmass_out > 10**6)])[0]) +'/'+str(np.shape(GW_tmass)[0]))
            file.write('\nFraction of (IMBH-SMBH)-IMBH tertiaries:         '+str(np.shape(GW_tmass[(GW_tmass_out < 10**6) & (GW_tmass_out > 10**6)])[0]) +'/'+str(np.shape(GW_tmass)[0]))
            file.write('\nFraction of (IMBH-IMBH)-IMBH tertiaries:         '+str(np.shape(GW_tmass[(GW_tmass_in < 10**6) & (GW_tmass_out < 10**6)])[0]) +'/'+str(np.shape(GW_tmass)[0]))

        return bin_formed, ter_formed, bsys_time, tsys_time, semi_major_bin, semi_major_ter
    
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
        
        dedt = [[ ], [ ]]
        dadt = [[ ], [ ]]
        for j in range(2):
            for sublist in self.dedt[j]:
                for item_ in sublist:
                    dedt[j].append(float(item_))
            for sublist in self.dadt[j]:
                for item_ in sublist:
                    dadt[j].append(item_)
            
            dedt[j] = np.asarray(dedt[j])
            dadt[j] = np.asarray(dadt[j])
            
        p21ymin = 6*np.log10(np.nanmin(dedt[0]))
        p21ymax = 6*np.log10(np.nanmax(dedt[0]))
        p21xmin = 6*np.log10(np.nanmin(dadt[0]))
        p21xmax = 6*np.log10(np.nanmax(dadt[0]))
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

        hardb_idx = [ ]
        softb_idx = [ ]
        hardt_idx = [ ]
        softt_idx = [ ]
        for int_ in range(1):
            bhidx = np.where(self.hard_bin[int_] > 0)[0]
            bsidx = np.where(self.hard_bin[int_] < 0)[0]
            thidx = np.where(self.hard_ter[int_] > 0)[0]
            tsidx = np.where(self.hard_ter[int_] < 0)[0]
            hardb_idx[int_].append(bhidx)
            softb_idx[int_].append(bsidx)
            hardt_idx[int_].append(thidx)
            softt_idx[int_].append(tsidx)

        ####### CALCULATE FOR TERTIARY, HARD BINARY, SOFT BINARY AND ALL + NOTE tGW < tH
        for int_ in range(1):
            int_ += 1
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


print('...sustainable_bintert_plotters...')
cst = sustainable_sys()
#cst.new_data_extractor()
cst.combine_data()
cst.system_formation_plotter()
#cst.GW_emissions()