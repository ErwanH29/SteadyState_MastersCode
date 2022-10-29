from amuse.lab import *
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
        
    def chaos_extract(self, dirC):
        """
        Function to extract data from simulations who end in ejection
        """
        self.chaos_ini_parti_data, self.chaos_fin_parti_data, self.chaos_number_mergers, self.chaos_cumulative_mm, self.chaos_simulated_end, \
        self.chaos_ejected_parti, self.chaos_stab_time_data, self.chaos_init_dist_data, self.chaos_init_mass_data, \
        self.chaos_inj_mass_data, self.chaos_eje_mass_data = stats_chaos_extractor(dirC)  

        return self.chaos_ini_parti_data, self.chaos_fin_parti_data, self.chaos_number_mergers, self.chaos_cumulative_mm, \
               self.chaos_simulated_end, self.chaos_ejected_parti, self.chaos_stab_time_data, self.chaos_init_dist_data, \
               self.chaos_init_mass_data, self.chaos_inj_mass_data, self.chaos_eje_mass_data

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
        ax.set_ylabel(r'$t_{eject}$ [Myr]')
        ax.set_title(r'Spread in Stability Time')
        
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
        y_labs = [r'$t_{eject}$ [Myr]', r'Percentage of Survival [%]']
        title = [r'Black Hole Population vs. Stability Time', r'Black Hole Population vs. Survivability Rate']
    
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

        print(avgt_arr)
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
            xints = [i for i in range(1+int(max(pops_f))) if i %10 == 0]
        if delim == 'B':
            xints = [i for i in range(1+int(max(pops_f)))]
            
        self.mean_plots([ax], pop_arr, xints, y_max)
        plot_ini.tickers(ax)
        ax.xaxis.labelpad = 30
        ax.legend()

    def stable_extract(self, dirS):
        """
        Function to extract data  from simulations who end with stable state
        """
        self.ini_parti_data, self.inj_event_data, self.merge_no_data, self.mergerm_data, self.simulated_end, \
        self.initial_dist, self.cluster_rad, self.init_parti_m, self.trelax_data = stats_stable_extractor(dirS)

        return self.ini_parti_data, self.inj_event_data, self.merge_no_data, self.mergerm_data, self.simulated_end, \
               self.initial_dist, self.cluster_rad, self.init_parti_m, self.trelax_data

    def massdep_plotter(self):
        """
        Function to plot stability time for constant distances
        """

        def log_fit(xval, slope, yint):
            return (slope) /( (xval)*np.log(xval))**1 + yint

        plot_ini = plotter_setup()

        dirH = 'data/Hermite/no_addition/chaotic_simulation/*'
        dirG = 'data/GRX/no_addition/chaotic_simulation/*'

        chaos_ini_parti_data, chaos_fin_parti_data, chaos_number_mergers, chaos_cumulative_mm, chaos_simulated_end, \
        chaos_ejected_parti, chaos_stab_time_data, chaos_init_dist_data, chaos_init_mass_data, \
        chaos_inj_mass_data, chaos_eje_mass_data = self.chaos_extract(dirH)

        chaos_ini_parti_data_GRX, chaos_fin_parti_data_GRX, chaos_number_mergers_GRX, chaos_cumulative_mm_GRX, chaos_simulated_end_GRX, \
        chaos_ejected_parti_GRX, chaos_stab_time_data_GRX, chaos_init_dist_data_GRX, chaos_init_mass_data_GRX, \
        chaos_inj_mass_data_GRX, chaos_eje_mass_data_GRX = self.chaos_extract(dirG)

        in_dist = np.unique(chaos_init_dist_data)
        in_mass = np.unique(chaos_init_mass_data, axis=0)

        for dist_ in in_dist:
            init_dist_idx_chaos = np.where((np.asarray(chaos_init_dist_data) == dist_))
            init_dist_idx_chaos_GRX = np.where((np.asarray(chaos_init_dist_data_GRX) == dist_))

            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111)

            tot_pop  = [ ]   
            iter = -1   
            coloursycler = cycle(colour_picker())
            for mass_ in in_mass:
                iter += 1
                pop_size, pop_samples, fin_parti, stab_time = self.index_extractor(chaos_init_mass_data, chaos_fin_parti_data, 
                                                                                   chaos_stab_time_data, init_dist_idx_chaos, mass_)
                pop_size_GRX, pop_samples_GRX, fin_parti_GRX, stab_time_GRX = self.index_extractor(chaos_init_mass_data_GRX, chaos_fin_parti_data_GRX, 
                                                                                   chaos_stab_time_data_GRX, init_dist_idx_chaos_GRX, mass_)

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

                full_simul = []
                for pop_, samp_ in zip(pop_size, pop_samples):
                    N_parti = np.argwhere(fin_parti == pop_)
                    N_parti_avg.append(np.mean(stab_time[N_parti]))
                    idx = np.where(stab_time[N_parti] == 100)[0]
                    ratio = len(idx)/len(stab_time[N_parti])
                    full_simul.append(ratio)
                N_parti_avg = np.asarray(N_parti_avg)
                print('For Hermite, # of full simulations per population: ', pop_size.flatten(), full_simul)

                full_simul = []
                for pop_, samp_ in zip(pop_size_GRX, pop_samples_GRX):
                    N_parti = np.argwhere(fin_parti_GRX == pop_)
                    N_parti_avg_GRX.append(np.mean(stab_time_GRX[N_parti]))
                    idx = np.where(stab_time_GRX[N_parti] == 100)[0]
                    ratio = len(idx)/len(stab_time_GRX[N_parti])
                    full_simul.append(ratio)
                N_parti_avg_GRX = np.asarray(N_parti_avg_GRX)
                print('For GRX, # of full simulations per population: ', pop_size_GRX.flatten(), full_simul)

                ax.scatter(pop_size, N_parti_avg, color = 'red', edgecolor = 'black', zorder = 3,
                           label = r'Hermite')                
                ax.scatter(pop_size_GRX, N_parti_avg_GRX, edgecolor = 'black', color = 'blue', 
                           zorder = 3, label = r'Hermite GRX')

            for j, xpos in enumerate(pop_size):
                ax.text(pop_size[j][0], -0.1*max(N_parti_avg), '# Ejec.\n'+'Hermite: '+str('{:.0f}'.format(pop_samples[j][0])), fontsize = 'xx-small', ha = 'center' )
            for j, xpos in enumerate(pop_size_GRX):
                ax.text(pop_size_GRX[j][0], -0.1*max(N_parti_avg_GRX), '# Ejec.\n'+'Hermite: '+str('{:.0f}'.format(pop_samples_GRX[j][0])), fontsize = 'xx-small', ha = 'center' )
            
            xints = [i for i in range(1+int(max(pop_size))) if i % 10 == 0]
            self.mean_plots([ax], pop_size, xints, max(N_parti_avg))
            plot_ini.tickers(ax)

            p0 = (2000, .1, 50)
            pop_size = np.array([float(i) for i in pop_size])
            N_parti_avg = np.array([ float(i) for i in N_parti_avg])

            p0 = (6,  0.50)
            params, cv = scipy.optimize.curve_fit(log_fit, pop_size, N_parti_avg, p0)
            slope, intercept = params
            red_slope = str('{:.2f}'.format(slope))
            xtemp = np.linspace(8, 105, 1000)
            ytemp = [log_fit(i, slope, intercept) for i in xtemp]
            ax.plot(xtemp, ytemp, zorder = 1, color = 'black', ls = '-.')
            ax.text(8, 1.5, r'$t_{{surv}} \approx \frac{{{}}}{{N\lnN}}$'.format(red_slope)+ ' Myr')

            ax.set_xlim(5,105)
            ax.legend()
            plt.savefig('figures/const_population_stab_time_equal_dist_'+str(dist_)+'_mean.pdf', dpi = 300, bbox_inches='tight')

        for dist_ in in_dist:
            for mass_ in in_mass: 
                y_max = []

                plt.clf()
                fig, ax = plt.subplots(figsize=(8,6))

                init_dist_idx = np.where((chaos_init_dist_data == dist_))
                init_mass = chaos_init_mass_data[init_dist_idx]
                fin_parti = chaos_fin_parti_data[init_dist_idx]
                stab_time = chaos_stab_time_data[init_dist_idx]

                tot_pop = [ ]
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
                    
                    ax.text(85, 0.96*(max(np.add(N_parti_avg, std))), r'$r_{SMBH}=$'+str(dist_)+' pc\n'+r'$m_{i} =$ '+str(mass_[0])+r' $M_\odot$')
                    xints = [i for i in range(1+int(max(tot_pop))) if i%10 == 0]
                    self.error_plots(ax, xints, tot_pop, N_parti_avg, std)
                    plot_ini.tickers(ax)

                    ax.set_xlim(5, 105)
                    plt.savefig('figures/const_pop_chaotic_stab_time_equal_dist_'+str(dist_)+'_err_mass_'+str(mass_)+'.pdf', dpi = 300, bbox_inches='tight')

cst = stability_plotters()
cst.massdep_plotter()


