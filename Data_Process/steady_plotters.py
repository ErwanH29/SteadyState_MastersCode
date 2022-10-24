from amuse.lab import *
from parti_initialiser import *
from file_logistics import *
from spatial_plotters import *
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import scipy.optimize
from itertools import cycle
from scipy import stats

class stability_plotters(object):
    """
    Class to extract data and setup plots for simulations where IMBH particles are added over time
    """
    def chaos_extract(self, dirC):
        """
        Function to extract data from simulations who end in ejection
        """
        self.chaos_ini_parti_data, self.chaos_fin_parti_data, self.chaos_number_mergers, self.chaos_cumulative_mm, self.chaos_simulated_end, \
        self.chaos_ejected_parti, self.chaos_stab_time_data, self.chaos_init_dist_data, self.chaos_cluster_radius, self.chaos_init_mass_data, \
        self.chaos_inj_mass_data, self.chaos_eje_mass_data, self.chaos_reltime_data = stats_chaos_extractor(dirC)  

        return self.chaos_ini_parti_data, self.chaos_fin_parti_data, self.chaos_number_mergers, self.chaos_cumulative_mm, \
               self.chaos_simulated_end, self.chaos_ejected_parti, self.chaos_stab_time_data, self.chaos_init_dist_data, \
               self.chaos_cluster_radius, self.chaos_init_mass_data, self.chaos_inj_mass_data, \
               self.chaos_eje_mass_data, self.chaos_reltime_data

    def stable_extract(self, dirS):
        """
        Function to extract data  from simulations who end with stable state
        """
        self.ini_parti_data, self.inj_event_data, self.merge_no_data, self.mergerm_data, self.simulated_end, \
        self.initial_dist, self.cluster_rad, self.init_parti_m, self.trelax_data = stats_stable_extractor(dirS)

        return self.ini_parti_data, self.inj_event_data, self.merge_no_data, self.mergerm_data, self.simulated_end, \
               self.initial_dist, self.cluster_rad, self.init_parti_m, self.trelax_data

    def axis2_filter(self, init_mass, iparti_data, inj_event_data, no_mergers, tot_pop, pop_samples, vals, indices):
        imass_surv = init_mass[indices]
        ini_parti = iparti_data[indices]
        inj_event = inj_event_data[indices]
        merge_event = no_mergers[indices]

        mass_arrays = np.where((imass_surv == vals).all(1))[0]  #Indices with the correct mass column
        ini_parti = ini_parti[mass_arrays]
        inj_event = inj_event[mass_arrays]
        merge_event = merge_event[mass_arrays]
        pop_sizeS, pop_samplesS = np.unique(ini_parti, return_counts=True)
        pop_sizeS = np.array([pop_+merge_-inj_ for pop_, merge_, inj_ in zip(pop_sizeS, merge_event, inj_event)])
        tot_popS = np.concatenate((pop_sizeS, tot_pop))
        tot_popS = np.unique(tot_popS)
        pop_samplesS = list(pop_samplesS)
        pop_sizeS = list(pop_sizeS)

        surv_rate = [ ]
        for i in range(len(tot_pop)):
            for j in range(len(pop_sizeS)):
                no_surv = True
                if tot_pop[i] == pop_sizeS[j]:
                    no_surv = False
                    surv_rate.append(100 * pop_samplesS[j]/(pop_samplesS[j]+pop_samples[i]))
                    break
            if (no_surv):
                surv_rate.append(100)
                pop_samplesS.insert(j, 1)
                pop_sizeS.insert(j, tot_pop[i])
        for j in range(len(pop_sizeS)):
            only_surv = True
            for i in range(len(tot_pop)):
                if pop_sizeS[j] == tot_pop[i]:
                    only_surv = False
                    break
            if (only_surv):
                surv_rate.insert(j, 0)
                pop_samplesS.insert(j, 1)
                pop_sizeS.insert(j, tot_pop[i])
        pop_samplesS = np.array(pop_samplesS)
        pop_sizeS = np.array(pop_sizeS)

        return pop_samplesS, pop_sizeS, tot_popS, surv_rate

    def error_plots(self, ax, xints, tot_pop, Navg, std):
        """
        Function to setup the error plots
        
        Input:
        ax:      The axis wanting to manipulate
        xints:   Array giving equal spacings based on range of sim. pop.
        tot_pop: Array holding of the range of simulated populations
        Navg:    Array with the average ejected times depending on population
        std:     The standard deviation for each simulation
        """

        ax.set_xlabel(r'Number of IMBH [$N$]')
        ax.set_xticks(xints)
        ax.set_xlim(2.5,max(tot_pop)+1)
        ax.set_ylim(0, 1.1*max(np.add(Navg, std)))
        ax.set_ylabel(r'Ejection Time [Myr]')
        ax.set_title(r'Chaotic Black Hole Population vs. Stability Time')
        
        return ax

    def index_extractor(self, init_mass, final_parti_data, stab_time_data, indices, vals):
        """
        Function to extract indices of depending on the wanted plot (function of mass/distance)
        
        Inputs:
        init_mass:         The set of initial masses
        final_parti_data:  The data of the final configuration
        stab_time_data:    The stability time
        indices:           The indices corresponding to the fashion you wish to separate the plots
        vals:              The variable which will compose different coloured lines in the plot
        """

        mass_array = init_mass[indices]
        final_part = final_parti_data[indices]
        stab_time = stab_time_data[indices]

        filtered_m = np.where((mass_array == vals).all(1))[0]
        filt_finparti = final_part[filtered_m]
        filt_stabtime = stab_time[filtered_m]
        pop_size, pop_samples = np.unique(filt_finparti, return_counts = True)

        return pop_size, pop_samples, final_part, stab_time

    def mean_plots(self, axis, pop, xints, ydata):
        """
        Function to set up the various stability time plotters
        
        Inputs:
        ax:     Array containing axes in plots
        pop:    The population, N, of the data
        xints:  The integer spacing for which to place the x-labels
        ydata:  The y-values for the plot
        """

        mtick_formatter = [mtick.FormatStrFormatter('%0.1f'), mtick.FormatStrFormatter('%0.0f')]
        ylims = [[0, 1.2*(ydata)], [0, 105]]
        y_labs = [r'Ejection Time [Myr]', r'Percentage of Survival [%]']
        title = [r'Black Hole Population vs. Average Ejection Time', r'Black Hole Population vs. Survivability Rate']
    
        for i in range(len(axis)):
            axis[i].set_xlabel(r'Number of IMBH [$N$]')
            axis[i].set_ylabel(y_labs[i])
            axis[i].set_xticks(xints)
            axis[i].set_xlim(2.5,max(pop)+1)
            axis[i].yaxis.set_major_formatter(mtick_formatter[i])
            axis[i].set_ylim(ylims[i])
            axis[i].set_title(title[i])
    
    def plot_splitter(self, pop_arr, samp_arr, avgt_arr, ax, delim):
        """
        Makes two separate plots depending on population
        
        Inputs:
        pop_arr:  The array with complete population samples
        samp_arr: The array denoting number of simulations conducted per pop.
        avgt_arr: The array corresponding to the average simulation time per pop.
        ax:       The axis to plot the figure on
        delim:    Delimiter splitting the plots ([A] for above/lower bound, [B] for upper bound)
        """

        plot_ini = plotter_setup()
        if delim == 'A':
            indx = np.argwhere(pop_arr > 9)
        if delim == 'B':
            indx = np.argwhere(pop_arr < 10)

        samps_f = samp_arr[indx[:,0]]
        pops_f = pop_arr[indx[:,0]]
        y_max = max(avgt_arr[indx[:,0]])
        for j, xpos in enumerate(pops_f):
            ax.text(pops_f[j][0], -0.15*(y_max), '# Ejec.\n'+'Hermite: '+str('{:.0f}'.format(samps_f[j][0])), fontsize = 'xx-small', ha = 'center' )
            #else:
            #   ax.text(xpos, -0.12*max(N_parti_avg)*(1+0.6*iter), 'Set '+str(iter)+': '+str(pop_samples[j][0]), fontsize = 'xx-small', ha = 'center' )
        if delim == 'A':
            xints = [i for i in range(1+int(max(pops_f))) if i %5 == 0]
        if delim == 'B':
            xints = [i for i in range(1+int(max(pops_f)))]
            
        self.mean_plots([ax], pop_arr, xints, y_max)
        plot_ini.tickers(ax)
        ax.xaxis.labelpad = 30
        ax.legend()

    def distdep_plotter(self, no_axis, int_string):
        """
        Function to plot the steady time with various lines corresponding to different distances.
        
        Inputs:
        no_axis: The number of axis wished to plot (1 = no added particles, 2 = particle flux)
        """
        
        plot_ini = plotter_setup()
        
        if no_axis == 1:
            dirH = 'data/Hermite/GC/no_addition/chaotic_simulation/*'
            dirG = 'data/GRX/GC/no_addition/chaotic_simulation/*'
        else:
            dirH = 'data/Hermite/GC/chaotic_simulation/*'
            dirG = 'data/GRX/GC/chaotic_simulation/*'

        chaos_ini_parti_data, chaos_fin_parti_data, chaos_number_mergers, chaos_cumulative_mm, chaos_simulated_end, \
        chaos_ejected_parti, chaos_stab_time_data, chaos_init_dist_data, chaos_cluster_radius, chaos_init_mass_data, \
        chaos_inj_mass_data, chaos_eje_mass_data, chaos_reltime_data = self.chaos_extract(dirH)

        chaos_ini_parti_data_GRX, chaos_fin_parti_data_GRX, chaos_number_mergers_GRX, chaos_cumulative_mm_GRX, chaos_simulated_end_GRX, \
        chaos_ejected_parti_GRX, chaos_stab_time_data_GRX, chaos_init_dist_data_GRX, chaos_cluster_radius_GRX, chaos_init_mass_data_GRX, \
        chaos_inj_mass_data_GRX, chaos_eje_mass_data_GRX, chaos_reltime_data_GRX = self.chaos_extract(dirG)

        if no_axis == 2:
            ini_parti_data, inj_event_data, merge_no_data, mergerm_data, simulated_end, \
            initial_dist, cluster_rad, init_parti_m, trelax_data = self.stable_extract('data/Hermite/GC/stable_simulation/*')

            ini_parti_data_GRX, inj_event_data_GRX, merge_no_data_GRX, mergerm_data_GRX, simulated_end_GRX, \
            initial_dist_GRX, cluster_rad_GRX, init_parti_m_GRX, trelax_data_GRX = self.stable_extract('data/GRX/GC/stable_simulation/*')

        in_dist = np.unique(chaos_init_dist_data)
        in_mass = np.unique(chaos_init_mass_data, axis=0)

        for mass_ in in_mass: #For every initial mass, we will plot a graph depending on rSMBH
            init_mass_idx = np.where((chaos_init_mass_data == mass_).all(1))[0] #Find indices where data files correspond to the correct initial masses
            init_mass_idx_GRX = np.where((chaos_init_mass_data_GRX == mass_).all(1))[0]

            if no_axis == 1:
                fig = plt.figure(figsize=(12.5, 8))
                ax11 = fig.add_subplot(221)
                ax12 = fig.add_subplot(222)

            if no_axis == 2:
                init_mass_idx_stab = np.where((init_parti_m == mass_).all(1))[0]
                init_mass_idx_stab_GRX = np.where((init_parti_m_GRX == mass_).all(1))[0]

                fig = plt.figure(figsize=(9, 14))
                ax = fig.add_subplot(211)
                ax2 = fig.add_subplot(212)
                tot_popS = [ ]

            tot_pop  = [ ]
            iter = -1
            colourcycler = cycle(colour_picker())
            for dist_ in in_dist:
                iter += 1
                fin_parti = chaos_fin_parti_data[init_mass_idx] #Filter the needed data based on the data files satisfying condition
                stab_time = chaos_stab_time_data[init_mass_idx]
                init_dist = chaos_init_dist_data[init_mass_idx]
                dist_arrays = np.where(init_dist == dist_) #Find indices corresponding to the correct distances. 
                fin_parti = fin_parti[dist_arrays] #Here we split the data relative to their distances.
                stab_time = stab_time[dist_arrays]
                pop_size, pop_samples = np.unique(fin_parti, return_counts=True) #Count the number of unique final populations

                fin_parti_GRX = chaos_fin_parti_data_GRX[init_mass_idx_GRX] 
                stab_time_GRX = chaos_stab_time_data_GRX[init_mass_idx_GRX]
                init_dist_GRX = chaos_init_dist_data_GRX[init_mass_idx_GRX]
                dist_arrays_GRX = np.where(init_dist_GRX == dist_)
                fin_parti_GRX = fin_parti[dist_arrays_GRX]         
                stab_time_GRX = stab_time[dist_arrays_GRX]
                pop_size_GRX, pop_samples_GRX = np.unique(fin_parti_GRX, return_counts=True)
               
                
                if no_axis == 2:
                    ini_parti = ini_parti_data[init_mass_idx_stab]
                    inj_event = inj_event_data[init_mass_idx_stab]
                    merge_event = merge_no_data[init_mass_idx_stab]
                    init_distS = initial_dist[init_mass_idx_stab]

                    dist_arraysS = np.where(init_distS == dist_)
                    ini_parti = ini_parti[dist_arraysS]
                    inj_event = inj_event[dist_arraysS]
                    merge_event = merge_event[dist_arraysS]
                    pop_sizeS, pop_samplesS = np.unique(ini_parti, return_counts=True)
                    pop_sizeS = np.array([pop_+merge_-inj_ for pop_, merge_, inj_ in zip(pop_sizeS, merge_event, inj_event)])
                    tot_popS = np.concatenate((pop_sizeS, pop_size))
                    tot_popS = np.unique(tot_popS)
                    pop_samplesS = list(pop_samplesS)
                    pop_sizeS = list(pop_sizeS)

                    surv_rate = [ ]
                    for i in range(len(pop_size)):
                        for j in range(len(pop_sizeS)):
                            no_surv = True
                            if pop_size[i] == pop_sizeS[j]:
                                no_surv = False
                                surv_rate.append(100 * pop_samplesS[j]/(pop_samplesS[j]+pop_samples[i]))
                                break
                        if (no_surv):
                            surv_rate.append(100)
                            pop_samplesS.insert(j, 1)
                            pop_sizeS.insert(j, pop_size[i])
                    for j in range(len(pop_sizeS)):
                        only_surv = True
                        for i in range(len(pop_size)):
                            if pop_sizeS[j] == pop_size[i]:
                                only_surv = False
                                break
                        if (only_surv):
                            surv_rate.insert(j, 0)
                            pop_samplesS.insert(j, 1)
                            pop_sizeS.insert(j, pop_size[i])
                    pop_samplesS = np.array(pop_samplesS)
                    pop_sizeS = np.array(pop_sizeS)

                N_parti_avg = [ ]
                pop_id = np.argwhere(pop_size > 2)  #Only look at samples with N>2
                pop_size = pop_size[pop_id]
                pop_samples = pop_samples[pop_id]

                N_parti_avg_GRX = [ ]
                pop_id_GRX = np.argwhere(pop_size_GRX > 2)
                pop_size_GRX = pop_size_GRX[pop_id_GRX]
                pop_samples_GRX = pop_samples_GRX[pop_id_GRX]
                
                colours = next(colourcycler)
                tot_pop.append(max(pop_size))
                if no_axis == 2:
                    pop_idS = np.argwhere(pop_sizeS > 2)
                    pop_sizeS = pop_sizeS[pop_idS]
                    pop_samplesS = pop_samplesS[pop_idS]

                    pop_idS_GRX = np.argwhere(pop_sizeS_GRX > 2)
                    pop_sizeS_GRX = pop_sizeS_GRX[pop_idS_GRX]
                    pop_samplesS_GRX = pop_samplesS_GRX[pop_idS_GRX]

                for pop_, samp_ in zip(pop_size, pop_samples):
                    N_parti = np.argwhere(fin_parti == pop_)
                    N_parti_avg.append(np.mean(stab_time[N_parti]))
                N_parti_avg = np.asarray(N_parti_avg)

                for pop_, samp_ in zip(pop_size_GRX, pop_samples_GRX):
                    N_parti = np.argwhere(fin_parti_GRX == pop_)
                    N_parti_avg_GRX.append(np.mean(stab_time_GRX[N_parti]))
                N_parti_avg_GRX = np.asarray(N_parti_avg_GRX)                     
                
                for ax_ in [ax11, ax12]:
                    ax_.scatter(pop_size, N_parti_avg, color = 'red', edgecolor = 'black', label = r'Hermite')
                    ax_.scatter(pop_size_GRX, N_parti_avg_GRX, color = 'blue', edgecolor = 'black', label = r'Hermit GRX')
                
                self.plot_splitter(pop_size, pop_samples, N_parti_avg, ax11, 'B')
                self.plot_splitter(pop_size, pop_samples, N_parti_avg, ax12, 'A')

                #self.plot_splitter(pop_size_GRX, pop_samples_GRX, N_parti_avg_GRX, ax11, 'B')
                #self.plot_splitter(pop_size_GRX, pop_samples_GRX, N_parti_avg_GRX, ax12, 'A')

            minorder = int('{:.0f}'.format(np.log10(mass_[0])))
            maxorder = int('{:.0f}'.format(np.log10(mass_[1])))

            if no_axis == 2:
                slope, intercept, r_value, p_value, std_err = stats.linregress(tot_popS, surv_rate)
                xtemp = np.linspace(2.6, max(tot_popS)+1)
                ytemp = [slope*xval_ + intercept for xval_ in xtemp]
                ax2.plot(xtemp, ytemp, color = 'black', zorder = 1)
                ax2.text(3.2, 10, r'Line of Best Fit: $S = $'+str('{:.2f}'.format(slope))+r'$N$')
                self.mean_plots([ax, ax2], tot_pop, xints, y_max)
                for ax_ in [ax, ax2]:
                    plot_ini.tickers(ax_)
                    ax_.xaxis.labelpad = 25
                ax.legend()
                plt.savefig('figures/chaotic_stab_time_equal_mass_'+str(mass_)+'_mean.pdf', dpi = 300, bbox_inches='tight')
            else:
                ax11.set_xlim(2.5,9.5)
                ax12.set_xlim(9.5,35.5)
                plt.savefig('figures/const_population_stab_time_equal_mass_'+str(mass_)+'_mean.pdf', dpi = 300, bbox_inches='tight')

        for mass_ in in_mass:      #For every initial mass and distance we will plot a separate graph. This time it includes std spread
            for dist_ in in_dist:  #For three different distances, every mass has 3 unique plots with errors.
                y_max = []
                
                plt.clf()
                fig, ax = plt.subplots()
                
                init_mass_idx = np.where((self.chaos_init_mass_data == mass_).all(1))[0]  #Make sure to use only data who has corresponding mass
                fin_parti = self.chaos_fin_parti_data[init_mass_idx]   #Filter data through indices satisfying the mass requirement
                stab_time = self.chaos_stab_time_data[init_mass_idx]
                init_dist = self.chaos_init_dist_data[init_mass_idx]
                dist_arrays = np.where(init_dist == dist_)  #Find indices corresponding to the correct distances. Here we split

                colours = next(colourcycler)
                N_parti_avg = [ ]
                std = [ ]
                
                fin_parti = fin_parti[dist_arrays]  #Filtering the filtered data with the correct distances
                stab_time = stab_time[dist_arrays]

                pop_size, pop_samples = np.unique(fin_parti, return_counts=True)  #Count the number of unique final populations
                
                if len(pop_size) == 0:
                    pass
                else:
                    pop_id = np.argwhere(pop_size > 2)
                    pop_size = pop_size[pop_id]
                    pop_samples = pop_samples[pop_id]

                    for pop_, samp_ in zip(pop_size, pop_samples):  #For the unique populations, compute their average stab time
                        N_parti = np.argwhere(fin_parti == pop_)    #and average spread in data.
                        N_parti_avg.append(np.mean(stab_time[N_parti]))
                        std.append(np.std(stab_time[N_parti]))

                    y_max.append(max(np.add(N_parti_avg, std)))
                    ax.errorbar(pop_size, N_parti_avg, color = 'black', yerr=std, fmt = 'o')
                    ax.scatter(pop_size, np.add(N_parti_avg, std), marker = '_', color = 'black')
                    ax.scatter(pop_size, np.subtract(N_parti_avg, std), marker = '_', color = 'black')
                    ax.scatter(pop_size, N_parti_avg, color = 'black')

                    for j, xpos in enumerate(pop_size):
                        ax.text(xpos, -0.115*max(np.add(N_parti_avg, std)), 'Ejected\n'+str(pop_samples[j][0]), fontsize = 'xx-small', ha = 'center' )
                    
                    order = int('{:.0f}'.format(np.log10(mass_)))
                    xints = [i for i in range(1+int(max(tot_pop)))]
                    ax.text(8, 0.96*(max(np.add(N_parti_avg, std))), r'$m_{IMBH}=$'+str('{:.0f}'.format(mass_/10**order))+r'$\times 10^{}$'.format(order)+r' $M_\odot$'+'\n'+r'$r_{SMBH}=$'+str(dist_)+' pc')
                         
                    self.error_plots(ax, xints, tot_pop, N_parti_avg, std)
                    plot_ini.tickers(ax)  
                    if no_axis == 1:
                        plt.savefig('figures/const_pop_chaotic_stab_time_equal_mass_'+str(mass_)+'_err_dist_'+str(dist_)+'.pdf', dpi = 300, bbox_inches='tight')
                    if no_axis == 2:
                        plt.savefig('figures/chaotic_stab_time_equal_mass_'+str(mass_)+'_err_dist_'+str(dist_)+'.pdf', dpi = 300, bbox_inches='tight')

    def massdep_plotter(self, no_axis):
        """
        Function to plot stability time for constant distances
        """

        def log_fit(xval, slope, yint):
            return (slope) /( (xval)*np.log(xval))**1 + yint

        plot_ini = plotter_setup()
        if no_axis == 1:
            dirH = 'data/Hermite/GC/no_addition/chaotic_simulation/*'
            dirG = 'data/GRX/GC/no_addition/chaotic_simulation/*'
        else:
            dirH = 'data/Hermite/GC/chaotic_simulation/*'
            dirG = 'data/GRX/GC/chaotic_simulation/*'

        chaos_ini_parti_data, chaos_fin_parti_data, chaos_number_mergers, chaos_cumulative_mm, chaos_simulated_end, \
        chaos_ejected_parti, chaos_stab_time_data, chaos_init_dist_data, chaos_cluster_radius, chaos_init_mass_data, \
        chaos_inj_mass_data, chaos_eje_mass_data, chaos_reltime_data = self.chaos_extract(dirH)

        chaos_ini_parti_data_GRX, chaos_fin_parti_data_GRX, chaos_number_mergers_GRX, chaos_cumulative_mm_GRX, chaos_simulated_end_GRX, \
        chaos_ejected_parti_GRX, chaos_stab_time_data_GRX, chaos_init_dist_data_GRX, chaos_cluster_radius_GRX, chaos_init_mass_data_GRX, \
        chaos_inj_mass_data_GRX, chaos_eje_mass_data_GRX, chaos_reltime_data_GRX = self.chaos_extract(dirG)
        #print(chaos_fin_parti_data, chaos_stab_time_data, chaos_init_dist_data, chaos_init_mass_data)
        if no_axis == 2:
            ini_parti_data, inj_event_data, merge_no_data, mergerm_data, simulated_end, \
            initial_dist, cluster_rad, init_parti_m, trelax_data = self.stable_extract('data/Hermite/GC/stable_simulation/*')

            ini_parti_data_GRX, inj_event_data_GRX, merge_no_data_GRX, mergerm_data_GRX, simulated_end_GRX, \
            initial_dist_GRX, cluster_rad_GRX, init_parti_m_GRX, trelax_data_GRX = self.stable_extract('data/GRX/GC/stable_simulation/*')

        in_dist = np.unique(chaos_init_dist_data)
        in_mass = np.unique(chaos_init_mass_data, axis=0)

        for dist_ in in_dist:
            init_dist_idx_chaos = np.where((np.asarray(chaos_init_dist_data) == dist_))
            init_dist_idx_chaos_GRX = np.where((np.asarray(chaos_init_dist_data_GRX) == dist_))
            if no_axis == 1:
                fig = plt.figure(figsize=(12.5, 8))
                ax11 = fig.add_subplot(221)
                ax12 = fig.add_subplot(222)

            if no_axis == 2:
                init_dist_idx_stab = np.where((initial_dist == dist_))
                init_dist_idx_stab_GRX = np.where((initial_dist_GRX == dist_))
                
                fig = plt.figure(figsize=(9, 14))
                ax = fig.add_subplot(211)
                ax2 = fig.add_subplot(212)

            tot_pop  = [ ]   
            iter = -1   
            coloursycler = cycle(colour_picker())
            for mass_ in in_mass:
                iter += 1
                pop_size, pop_samples, fin_parti, stab_time = self.index_extractor(chaos_init_mass_data, chaos_fin_parti_data, 
                                                                                   chaos_stab_time_data, init_dist_idx_chaos, mass_)
                pop_size_GRX, pop_samples_GRX, fin_parti_GRX, stab_time_GRX = self.index_extractor(chaos_init_mass_data_GRX, chaos_fin_parti_data_GRX, 
                                                                                   chaos_stab_time_data_GRX, init_dist_idx_chaos_GRX, mass_)
                if no_axis == 2:
                    pop_samplesS, pop_sizeS, tot_popS, surv_rate = self.axis2_filter(init_parti_m, ini_parti_data, inj_event_data, 
                                                                                     merge_no_data, pop_size, pop_samples, mass_, 
                                                                                     init_dist_idx_stab)

                    pop_samplesS_GRX, pop_sizeS_GRX, tot_popS_GRX, surv_rate_GRX = self.axis2_filter(init_parti_m_GRX, ini_parti_data_GRX, inj_event_data_GRX, 
                                                                                                     merge_no_data_GRX, pop_size_GRX, pop_samples_GRX, mass_, 
                                                                                                     init_dist_idx_stab_GRX)
                N_parti_avg = [ ]
                pop_id = np.argwhere(pop_size > 2)
                pop_size = pop_size[pop_id]
                pop_samples = pop_samples[pop_id]

                N_parti_avg_GRX = [ ]
                pop_id_GRX = np.argwhere(pop_size_GRX > 2)
                pop_size_GRX = pop_size_GRX[pop_id_GRX]
                pop_samples_GRX = pop_samples_GRX[pop_id_GRX]

                colours = next(coloursycler)
                tot_pop.append(max(pop_size))
                if no_axis == 2:
                    pop_idS = np.argwhere(pop_sizeS > 2)
                    pop_sizeS = pop_sizeS[pop_idS]
                    pop_samplesS = pop_samplesS[pop_idS]

                    pop_idS_GRX = np.argwhere(pop_sizeS_GRX > 2)
                    pop_sizeS_GRX = pop_sizeS_GRX[pop_idS_GRX]
                    pop_samplesS_GRX = pop_samplesS_GRX[pop_idS_GRX]

                for pop_, samp_ in zip(pop_size, pop_samples):
                    N_parti = np.argwhere(fin_parti == pop_)
                    N_parti_avg.append(np.mean(stab_time[N_parti]))
                N_parti_avg = np.asarray(N_parti_avg)

                for pop_, samp_ in zip(pop_size_GRX, pop_samples_GRX):
                    N_parti = np.argwhere(fin_parti_GRX == pop_)
                    N_parti_avg_GRX.append(np.mean(stab_time_GRX[N_parti]))
                N_parti_avg_GRX = np.asarray(N_parti_avg_GRX)

                for ax_ in [ax11, ax12]:
                    ax_.scatter(pop_size, N_parti_avg, color = 'red', edgecolor = 'black', zorder = 3,
                                label = r'Hermite')                
                    ax_.scatter(pop_size_GRX, N_parti_avg_GRX, edgecolor = 'black', color = 'blue', 
                                zorder = 3, label = r'Hermite GRX')

                self.plot_splitter(pop_size, pop_samples, N_parti_avg, ax11, 'B')
                self.plot_splitter(pop_size, pop_samples, N_parti_avg, ax12, 'A')

                #self.plot_splitter(pop_size_GRX, pop_samples_GRX, N_parti_avg_GRX, ax11, 'B')
                #self.plot_splitter(pop_size_GRX, pop_samples_GRX, N_parti_avg_GRX, ax12, 'A')

                if no_axis == 2:
                    ax2.scatter(tot_popS, surv_rate, color = colours, edgecolors = 'black', zorder = 2)
                    for j, xpos in enumerate(pop_sizeS):
                        if iter == 0:
                            ax2.text(xpos, -0.12*max(surv_rate), '# Surv.\n'+'Set '+str(iter)+': '+str('{:.0f}'.format(pop_samplesS[j][0])), fontsize = 'xx-small', ha = 'center' )
                        else:
                            ax2.text(xpos, -0.12*max(surv_rate)*(1+0.6*iter), 'Set '+str(iter)+': '+str(pop_samplesS[j][0]), fontsize = 'xx-small', ha = 'center' )
                            
                    ax2.scatter(tot_popS_GRX, surv_rate_GRX, color = 'blue', edgecolors = 'black', zorder = 2)
                    for j, xpos in enumerate(pop_sizeS_GRX):
                        ax2.text(xpos, -0.15*max(surv_rate)*(1+1.2*iter), 'Set '+str(iter)+': '+str(pop_samplesS_GRX[j][0]), fontsize = 'xx-small', ha = 'center' )

            gc_code =globular_cluster()
            p0 = (2000, .1, 50)
            pop_size = np.array([float(i) for i in pop_size])
            N_parti_avg = np.array([ float(i) for i in N_parti_avg])

            p0 = (6,  0.50)
            params, cv = scipy.optimize.curve_fit(log_fit, pop_size, N_parti_avg, p0)
            slope, intercept = params
            red_slope = str('{:.2f}'.format(slope))
            xtemp = np.linspace(2.5, 50, 1000)
            ytemp = [log_fit(i, slope, intercept) for i in xtemp]
            for ax_ in [ax11, ax12]:
                ax_.plot(xtemp, ytemp, zorder = 1, color = 'black', ls = '-.')
            ax11.text(2.8, 0.5, r'$t_{{surv}} \approx \frac{{{}}}{{N\lnN}}$'.format(red_slope)+ ' Myr')
            
            if no_axis == 2:
                slope, intercept, r_value, p_value, std_err = stats.linregress(tot_popS, surv_rate)
                xtemp = np.linspace(2.6, max(tot_popS)+1)
                ytemp = [slope*xval_ + intercept for xval_ in xtemp]
                ax2.plot(xtemp, ytemp, color = 'black', zorder = 2)

                slope_GRX, intercept_GRX, r_value, p_value, std_err = stats.linregress(tot_popS_GRX, surv_rate_GRX)
                xtemp_GRX = np.linspace(2.6, max(tot_popS_GRX)+1)
                ytemp_GRX = [slope_GRX*xval_ + intercept_GRX for xval_ in xtemp_GRX]
                ax2.plot(xtemp_GRX, ytemp_GRX, color = 'black', zorder = 2)

                ax2.text(3.2, 10, r'Line of Best Fit: $S = $'+str('{:.2f}'.format(slope))+r'$N$')
                ax2.xaxis.labelpad = 25
                self.mean_plots([ax, ax2], tot_pop, xints, y_max)                
                for ax_ in [ax, ax2]:
                    plot_ini.tickers(ax_)
                ax.legend()
                plt.savefig('figures/chaotic_stab_time_equal_dist_'+str(dist_)+'_mean.pdf', dpi = 300, bbox_inches='tight')
            else:
                ax11.set_xlim(2.5,9.5)
                ax12.set_xlim(9.5,35.5)
                plt.savefig('figures/const_population_stab_time_equal_dist_'+str(dist_)+'_mean.pdf', dpi = 300, bbox_inches='tight')

        for dist_ in in_dist:
            for mass_ in in_mass: 
                y_max = []

                plt.clf()
                fig, ax = plt.subplots(figsize=(9,6))

                init_dist_idx = np.where((chaos_init_dist_data == dist_))
                init_mass = chaos_init_mass_data[init_dist_idx]
                fin_parti = chaos_fin_parti_data[init_dist_idx]
                stab_time = chaos_stab_time_data[init_dist_idx]

                tot_pop = [ ]
                colours = next(coloursycler)
                N_parti_avg = [ ]
                std = [ ]
                mass_arrays = np.where((init_mass == mass_).all(1))[0]  #Indices with the correct mass column

                fin_parti = fin_parti[mass_arrays]
                stab_time = stab_time[mass_arrays]
                pop_size, pop_samples = np.unique(fin_parti, return_counts=True)
                pop_id = np.argwhere(pop_size > 2)
                pop_size = pop_size[pop_id]

                if len(pop_size) == 0:
                    pass
                else:
                    pop_samples = pop_samples[pop_id]
                    tot_pop.append(max(pop_size))

                    for pop_, samp_ in zip(pop_size, pop_samples):
                        N_parti = np.argwhere(fin_parti == pop_)
                        N_parti_avg.append(np.mean(stab_time[N_parti]))
                        std.append(np.std(stab_time[N_parti]))

                    y_max.append(max(np.add(N_parti_avg, std)))
                    ax.errorbar(pop_size, N_parti_avg, color = 'black', yerr=std, fmt = 'o')
                    ax.scatter(pop_size, np.add(N_parti_avg, std), marker = '_', color = 'black')
                    ax.scatter(pop_size, np.subtract(N_parti_avg, std), marker = '_', color = 'black')
                    ax.scatter(pop_size, N_parti_avg, color = 'black')

                    #for j, xpos in enumerate(pop_size):
                    #    ax.text(xpos, -0.12*max(np.add(N_parti_avg, std)), '# Sim:\n'+str('{:.0f}'.format(pop_samples[j][0])), fontsize = 'xx-small', ha = 'center' )                   
                    ax.text(9, 0.96*(max(np.add(N_parti_avg, std))), r'$r_{SMBH}=$'+str(dist_)+' pc\n'+r'$m_{i} =$ '+str(mass_[0])+r' $M_\odot$')
                    xints = [i for i in range(1+int(max(tot_pop)))]
                    
                    self.error_plots(ax, xints, tot_pop, N_parti_avg, std)
                    plot_ini.tickers(ax)
                    ax.set_xlim(2.5,10.5)
                    if no_axis == 1:
                        plt.savefig('figures/const_pop_chaotic_stab_time_equal_dist_'+str(dist_)+'_err_mass_'+str(mass_)+'.pdf', dpi = 300, bbox_inches='tight')
                    if no_axis == 2:
                        plt.savefig('figures/chaotic_stab_time_equal_dist_'+str(dist_)+'_err_mass_'+str(mass_)+'.pdf', dpi = 300, bbox_inches='tight')


string = 'GRX'

#cst = stability_plotters()
#cst.massdep_plotter(1)
#cst.distdep_plotter(1, 'Hermite')
#spatial_plotter('Hermite')

#vej_plot = vejec_mass()
#vej_plot.vejec_syspop()
#vej_plot.vejec_sysmass()


