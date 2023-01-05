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

        
        dir = os.path.join('/home/erwanh/Desktop/SteadyStateBH/Data_Process/figures/steady_time/Sim_summary.txt')
        with open(dir) as f:
            line = f.readlines()

            popG = line[14][48:-2] 
            avgG = line[21][48:-2]
            avgG2 = line[22][3:-2] 
            popG_data = popG.split()
            avgG_data = avgG.split()
            avgG2_data = avgG2.split()
            popG = np.asarray([float(i) for i in popG_data])
            avgG = np.asarray([float(i) for i in avgG_data])
            avgG2 = np.asarray([float(i) for i in avgG2_data])
            avgG = np.concatenate((avgG, avgG2))

        for int_ in range(2):
            for file_ in range(len(filename[int_])):
                with open(filename[int_][file_], 'rb') as input_file:
                    print('Reading file', file_, ':', input_file)
                    file_size = os.path.getsize(filename[int_][file_])
                    if file_size < 2.8e9:
                        count += 1
                        data = pkl.load(input_file)
                        pop = 5*round(0.2*np.shape(data)[0])
                        idx = np.where(popG == pop)
                        col_len = int(min(np.round((avgG[idx])*10**3), np.shape(data)[1])-1)

                        if pop > 5 and pop <= 40:
                            for parti_ in range(np.shape(data)[0]):
                                integrator = ints[int_]

                                pop_bin = [ ]
                                pop_ter = [ ]

                                semi_NN_avg = [ ]
                                semi_NN_min = [ ]
                                semi_t_avg = [ ]
                                semi_t_min = [ ]

                                bsys_time = [ ]
                                bform_time = [ ]
                                bin_key = [ ]
                                new_bin_sys = [ ]
                                tsys_time = [ ]
                                tform_time = [ ]
                                ter_key = [ ]
                                new_ter_sys = [ ]

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
                                    temp_dedt.append((data.iloc[parti_][col_len][8][0] - data.iloc[parti_][2][8][0])/(np.shape(data)[1]-3))
                                    temp_dadt.append((data.iloc[parti_][col_len][7][0]**-1 - data.iloc[parti_][2][7][0]**-1).value_in((units.pc)**-1)/(np.shape(data)[1]-3))
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
                                                    bin_key.append(data.iloc[parti_][col_][6][1])
                                                    if data.iloc[parti_][col_][6][1] != bin_key[-1] or len(new_bin_sys) < 1:
                                                        new_bin_sys.append(1)
                                                    else:
                                                        new_bin_sys.append(-5)

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
                                                            ter_key.append([i for i in data.iloc[parti_][col_][6][2]][0])
                                                            if data.iloc[parti_][col_][6][2] != ter_key[-1] or len(new_ter_sys) < 2:
                                                                new_ter_sys.append(1)
                                                            else:
                                                                new_ter_sys.append(-5)

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

                                    bin_sys = np.shape(np.unique(bin_key))[0]
                                    ter_sys = np.shape(np.unique(ter_key))[0]

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
                                                
                                    bsys_time.append((len(np.unique(bsys_time)))/(col_))
                                    tsys_time.append((len(np.unique(tsys_time)))/(col_))
                                    dedt.append(np.mean(temp_dedt))
                                    dadt.append(np.mean(temp_dadt))

                                    path = '/home/erwanh/Desktop/SteadyStateBH/Data_Process/data/bin_hier_systems/'
                                    stab_tracker = pd.DataFrame()
                                    df_stabtime = pd.Series({'Integrator': integrator,
                                                            'Population': pop,
                                                            'Binary Pop.': pop_bin,
                                                            '# Binary Sys.': bin_sys,
                                                            'First bin. form': bform_time,
                                                            'Binary System Delineate': new_bin_sys,
                                                            'Tertiary Pop.': pop_ter,
                                                            '# Tertiary Sys.': ter_sys,
                                                            'First ter. form': tform_time,
                                                            'Tertiary System Delineate': new_ter_sys,
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
                                                            'Ter. Semi-major Min': semi_t_min,
                                                            'Total sim. length': col_ * 1000,
                                                            })
                                    stab_tracker = stab_tracker.append(df_stabtime, ignore_index = True)
                                    stab_tracker.to_pickle(os.path.join(path, 'IMBH_'+str(ints[int_])+'_system_data_indiv_parti_'+str(count)+'_'+str(parti_)+'_local2.pkl'))

    def array_rewrite(self, arr, arr_type, filt):
        """
        Function to rewrite array to manipulatable float format

        Inputs:
        arr:      The original array with the data
        new_arr:  The new, manipulatable, array
        arr_type: String stating whether array is nested or not 
        filt:     Boolean to filter out unwanted values
        """

        new_arr = [ ]
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

        GW_calcs = gw_calcs()

        self.integrator = [[ ], [ ]]
        self.pop = [[ ], [ ]]
        self.pop_bin = [[ ], [ ]]
        self.bin_sys = [[ ], [ ]]
        self.bform_time = [[ ], [ ]]
        self.new_bin_sys = [[ ], [ ]]
        self.pop_ter = [[ ], [ ]]
        self.ter_sys = [[ ], [ ]]
        self.tform_time = [[ ], [ ]]
        self.new_ter_sys = [[ ], [ ]]
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
        self.tot_sim = [[ ], [ ]]

        system_data = natsort.natsorted(glob.glob('data/bin_hier_systems/*'))
        for file_ in range(len(system_data)):
            with open(system_data[file_], 'rb') as input_file:
                data_file = pkl.load(input_file)
                if data_file.iloc[0][0] == 'Hermite':
                    int_ = 0
                else:
                    int_ = 1
                self.integrator[int_].append(data_file.iloc[0][0][0])
                self.pop[int_].append(int(data_file.iloc[0][1]))
                self.pop_bin[int_].append(data_file.iloc[0][2])
                self.bin_sys[int_].append(data_file.iloc[0][3])
                self.bform_time[int_].append(data_file.iloc[0][4])
                self.new_bin_sys[int_].append(data_file.iloc[0][5])
                self.pop_ter[int_].append(data_file.iloc[0][6])
                self.ter_sys[int_].append(data_file.iloc[0][7])
                self.tform_time[int_].append(data_file.iloc[0][8])
                self.new_ter_sys[int_].append(data_file.iloc[0][9])
                self.hard_bin[int_].append(data_file.iloc[0][10])
                self.hard_ter[int_].append(data_file.iloc[0][11])
                self.GW_freqbin[int_].append(data_file.iloc[0][12])
                self.GW_strainbin[int_].append(data_file.iloc[0][13])
                self.GW_timeb[int_].append(data_file.iloc[0][14])
                self.GW_bmass[int_].append(data_file.iloc[0][15])
                self.GW_freqter[int_].append(data_file.iloc[0][16])
                self.GW_strainter[int_].append(data_file.iloc[0][17])
                self.GW_timet[int_].append(data_file.iloc[0][18])
                self.GW_tmass[int_].append(data_file.iloc[0][19])
                self.dedt[int_].append(data_file.iloc[0][20])
                self.dadt[int_].append(data_file.iloc[0][21])
                self.semi_NN_avg[int_].append(data_file.iloc[0][22])
                self.semi_NN_min[int_].append(data_file.iloc[0][23])
                self.semi_t_avg[int_].append(data_file.iloc[0][24])
                self.semi_t_min[int_].append(data_file.iloc[0][25])
                self.tot_sim[int_].append(data_file.iloc[0][26])
                
        self.binary_systems = [[ ], [ ]]
        self.binary_occupation = [[ ], [ ]]
        self.tertiary_systems = [[ ], [ ]]
        self.tertiary_occupation = [[ ], [ ]]
        with open('figures/binary_hierarchical/output/system_summary.txt', 'w') as file:
            integrator = ['Hermite', 'GRX']
            for int_ in range(2):
                bform_stat = [ ]
                tform_stat = [ ]
                minGWtime = [ ]
                semi_NN_astat = [ ]
                semi_NN_mstat = [ ]
                semi_t_astat = [ ]
                semi_t_mstat = [ ]
                binHubble = [ ]
                totalbin = [ ]

                pop_arr = np.unique(self.pop[int_])
                file.write('Data for '+str(integrator[int_]+' in pc'))
                for pop_ in pop_arr:
                    idx = np.argwhere(self.pop[int_] == pop_).flatten()
                    bform_time = [ ]
                    tform_time = [ ]
                    GWb_time = [ ]
                    semi_NN_avg = [ ]
                    semi_NN_min = [ ]
                    semi_t_avg = [ ]
                    semi_t_min = [ ]
                    bin_occ = [ ]
                    ter_occ = [ ]

                    bin_data = 0
                    bin_data_merge = 0
                    for data_ in idx:
                        for item_ in self.bform_time[int_][data_]:
                            if item_ > 0:
                                bform_time.append(item_)
                        for item_ in self.tform_time[int_][data_]:
                            if item_ > 0:
                                tform_time.append(item_)
                        idx_nbin = np.argwhere(np.asarray(self.new_bin_sys[int_][data_]) > 0).flatten()
                        GW_btime_temp = np.asarray(self.GW_timeb[int_][data_])[idx_nbin-1]
                        bin_data += 1
                        if GW_btime_temp < (GW_calcs.tH).value_in(units.Myr):
                            bin_data_merge += 1
                            GWb_time.append(float(GW_btime_temp))

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
                        bin_occ.append(1000*(len(self.pop_bin[int_][data_]) / self.tot_sim[int_][data_]))
                        ter_occ.append(1000*(len(self.pop_ter[int_][data_]) / self.tot_sim[int_][data_]))


                    binHubble.append(bin_data_merge)
                    totalbin.append(bin_data)
                    
                    self.binary_systems[int_].append(bin_data/len(idx))
                    self.binary_occupation[int_].append(np.mean(bin_occ))
                    self.tertiary_systems[int_].append(len(np.asarray(self.semi_t_min[int_][data_])[np.asarray(self.semi_t_min[int_][data_]) > 0])/len(idx))
                    self.tertiary_occupation[int_].append(np.mean(ter_occ))

                    bform_stat.append('{:.7f}'.format(np.mean(bform_time)))
                    tform_stat.append('{:.7f}'.format(np.mean(tform_time)))
                    minGWtime.append(np.sort(GWb_time))
                    semi_NN_astat.append('{:.7f}'.format(np.mean(semi_NN_avg)))
                    semi_NN_mstat.append('{:.7f}'.format(np.mean(semi_NN_min)))
                    semi_t_astat.append('{:.7f}'.format(np.mean(semi_t_avg)))
                    semi_t_mstat.append('{:.7f}'.format(np.mean(semi_t_min)))

                file.write('\nAverage binary formation time [yrs]:             '+str(pop_arr)+' : '+str(bform_stat))
                file.write('\nAverage tertiary formation time [yrs]:           '+str(pop_arr)+' : '+str(tform_stat))
                file.write('\nFinal GW timescale (only tGW < tH) [Myr]:        '+str(pop_arr)+' : '+str(minGWtime))
                file.write('\nFraction of mergers within Hubble time:          '+str(pop_arr)+' : '+str(binHubble)+' / '+str(totalbin))
                file.write('\nAverage binary semi-major axis per population:   '+str(pop_arr)+' : '+str(semi_NN_astat))
                file.write('\nMinimum binary semi-major axis per population:   '+str(pop_arr)+' : '+str(semi_NN_mstat))
                file.write('\nAverage tertiary semi-major axis per population: '+str(pop_arr)+' : '+str(semi_t_astat))
                file.write('\nMinimum tertiary semi-major axis per population: '+str(pop_arr)+' : '+str(semi_t_mstat)+'\n\n')
                
        """with open('figures/binary_hierarchical/output/'+str(integrator[int_])+'bin_ter_systems.txt', 'w') as file:
            file.write(str(integrator[int_])+' BINARY DATA')
            file.write('\n# Binary systems at initialisation:             '+str(init_bin)+' / '+str(bin))
            file.write('\nFraction of IMBH-IMBH binaries:                 '+str(np.shape(GW_bmass[GW_bmass < 10**6])[0]) +'/'+str(np.shape(GW_bmass)[0]))
            file.write('\nFraction of hard binaries:                      '+str(np.sum(bhard[bhard > 0])) +'/'+str(np.shape(bhard)[0]))

            file.write('\n\nTERTIARY DATA')
            for pop_ in range(len(tform_time)):
                file.write('\nPopulation: '+str(ini_pop[pop_])+':         '+str(tform_time[pop_]/10**3)+' kyr')
            file.write('\n# Tertiary systems at initialisation:           '+str(init_ter)+' / '+str(ter))
            file.write('\nFraction of hard tertriary:                     '+str(np.sum(thard[thard > 0])) +'/'+str(np.shape(thard)[0]))

            file.write('\nFraction of mergers within Hubble time:         '+str(np.shape(GWt_mergertime[GWt_mergertime < tH.value_in(units.Myr)])[0]) +'/'+str(np.shape(GWt_mergertime)[0]))
            
            file.write('\nFraction of (IMBH-IMBH)-SMBH tertiaries:         '+str(np.shape(GW_tmass[(GW_tmass_in < 10**6) & (GW_tmass_out > 10**6)])[0]) +'/'+str(np.shape(GW_tmass)[0]))
            file.write('\nFraction of (IMBH-SMBH)-IMBH tertiaries:         '+str(np.shape(GW_tmass[(GW_tmass_out < 10**6) & (GW_tmass_out > 10**6)])[0]) +'/'+str(np.shape(GW_tmass)[0]))
            file.write('\nFraction of (IMBH-IMBH)-IMBH tertiaries:         '+str(np.shape(GW_tmass[(GW_tmass_in < 10**6) & (GW_tmass_out < 10**6)])[0]) +'/'+str(np.shape(GW_tmass)[0]))"""
    
    def system_formation_plotter(self):
        """
        Function to plot various 'sustainable system' plots
        """

        plot_ini = plotter_setup()
        mtick_formatter_y = mtick.FormatStrFormatter('%0.4f')
        mtick_formatter_x = mtick.FormatStrFormatter('%0.2f')
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
            
        p21ymin = 1.1*(np.nanmin((dedt[0])))
        p21ymax = 1.1*(np.nanmax((dedt[0])))
        p21xmin = 1.1*(np.nanmin((dadt[0])))
        p21xmax = 1.1*(np.nanmax((dadt[0])))
        normalise_p1 = plt.Normalize(0, max(max(self.binary_systems[0]), max(self.binary_systems[1])))
        normalise_p2 = plt.Normalize(10, 40)
    
        fig = plt.figure(figsize=(11, 6))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax_ = [ax1, ax2]
        for int_ in range(2):
            ini_pop = np.unique(self.pop[int_])
            ax_[int_].set_title(integrator[int_])
            ax_[int_].set_xlabel(r'IMBH Population [$N$]')
            ax_[int_].set_ylabel(r'$\log_{10}(t_{\rm{sys}} / t_{\rm{sim}})$')
            ax_[int_].set_ylim(-4, 0)
            colour_axes = ax_[int_].scatter(ini_pop, np.log10(self.binary_occupation[int_]), edgecolors  = 'black', c = (self.binary_systems[int_]), norm = (normalise_p1), label = 'Stable Binary')
            ax_[int_].scatter(ini_pop, np.log10(self.tertiary_occupation[int_]), edgecolors  = 'black', c = (self.tertiary_systems[int_]), norm = (normalise_p1), marker = 's', label = 'Stable Triple')
        plot_ini.tickers_pop(ax1, self.pop[0], 'Hermite')
        plot_ini.tickers_pop(ax2, self.pop[1], 'GRX')

        x_arr = np.linspace(5,110)
        y_arr1 = [0.02*i-3 for i in x_arr]
        y_arr2 = [0.05*i-3 for i in x_arr]
        ax1.plot(x_arr, y_arr1, color = 'black', linestyle = ':')
        ax2.plot(x_arr, y_arr2, color = 'black', linestyle = ':')
        ax2.legend()
        plt.colorbar(colour_axes, ax=ax2, label = r'$\langle N_{\rm{sys}} \rangle$ ')
        plt.savefig('figures/binary_hierarchical/sys_formation_N_plot.pdf', dpi=300, bbox_inches='tight')

        fig = plt.figure(figsize=(10, 6))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax = [ax1, ax2]
        ax1.set_ylabel(r'$\langle \dot{(1-e)} \rangle_{\rm{SMBH}}$ [Gyr$^{-1}$]')
        for ax_ in ax:
            plot_ini.tickers(ax_, 'plot')
            ax_.set_xlabel(r'$\langle \dot{a}^{-1} \rangle_{\rm{SMBH}}$ [pc$^{-1}$Gyr$^{-1}$]')            
            ax_.set_xlim(0.9*p21xmin, 1.1*p21xmax)
            ax_.set_ylim(0.9*p21ymin, 1.1*p21ymax)
            ax_.yaxis.set_major_formatter(mtick_formatter_y)        
            ax_.xaxis.set_major_formatter(mtick_formatter_x)
        for int_ in range(2):
            print(dadt, dedt)
            ax[int_].axvline(0, color = 'black', linestyle = ':')
            ax[int_].text(0.001, 0.0025, 'Hardening', rotation = 270)
            ax[int_].text(-0.004, 0.0025, 'Softening', rotation = 90)
            colour_axes = ax[int_].scatter((dadt[int_]), (dedt[int_]), 
                                           norm = normalise_p2, edgecolors = 'black', c = self.pop[int_])
            ax[int_].set_title(integrator[int_])

        plt.colorbar(colour_axes, ax = ax2, label = r'$N_{\rm{IMBH}}$')
        plt.savefig('figures/binary_hierarchical/simavg_dadt_dedt_plot.pdf', dpi=300, bbox_inches='tight')
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
cst.new_data_extractor()
cst.combine_data()
cst.system_formation_plotter()
#cst.GW_emissions()