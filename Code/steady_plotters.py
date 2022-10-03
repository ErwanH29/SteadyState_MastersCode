from sqlite3 import adapt
from amuse.lab import *
from parti_initialiser import *
from file_logistics import *
from spatial_plotters import *
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
from itertools import cycle
from scipy import stats

class stability_time(object):
    def __init__(self, dirC = 'data/chaotic_simulation/*', dirS = 'data/stable_simulation/*'):
        self.chaos_ini_parti_data, self.chaos_fin_parti_data, self.chaos_number_mergers, self.chaos_cumulative_mm, self.chaos_simulated_end, self.chaos_ejected_parti, \
        self.chaos_stab_time_data, self.chaos_init_dist_data, self.chaos_cluster_radius, self.chaos_init_mass_data, self.chaos_inj_mass_data, \
        self.chaos_eje_mass_data, self.chaos_reltime_data = stats_chaos_extractor(dirC)  

        self.ini_parti_data, self.inj_event_data, self.merge_no_data, self.mergerm_data, self.simulated_end, \
        self.initial_dist, self.cluster_rad, self.init_parti_m, self.trelax_data = stats_stable_extractor(dirS)

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
        ax.xaxis.labelpad = 20
        ax.set_xlabel(r'Number of IMBH [$N$]')
        ax.set_xticks(xints)
        ax.set_xlim(2.5,max(tot_pop)+1)
        ax.set_ylim(0, 1.1*max(np.add(Navg, std)))
        ax.set_ylabel(r'Ejection Time [Myr]')
        ax.set_title(r'Chaotic Black Hole Population vs. Stability Time')
        ax.xaxis.labelpad = 25
        return ax

    def stab_plotter_logistics(self, ax1, ax2, pop, xints, ydata):
        """
        Function to set up the various stability time plotters
        
        Inputs:
        ax:     The axis
        pop:    The population, N, of the data
        xints:  The integer spacing for which to place the x-labels
        ydata:  The y-values for the plot
        """

        for ax_ in [ax1, ax2]:
            ax_.set_xlabel(r'Number of IMBH [$N$]')
            ax_.set_xticks(xints)
            ax_.set_xlim(2.5,max(pop)+1)
        ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%0.3f'))
        ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%0.0f'))
        ax1.set_ylim(0, 1.2*max(ydata))
        ax2.set_ylim(0, 105)
        ax1.set_ylabel(r'Ejection Time [Myr]')
        ax2.set_ylabel(r'Percentage of Survival [%]')
        ax1.set_title(r'Chaotic Black Hole Population vs. Stability Time')
        ax2.set_title(r'Stable Black Hole Population Survivability Rate')

    def steadytime_distdep_plotter(self):
        """
        Function to plot the steady time with various lines corresponding to different
        SMBH Distances:
        """

        plot_ini = plotter_setup()

        in_mass = np.unique(self.chaos_init_mass_data[:,0], axis = 0) #Get the unique initial mass and distances array
        in_dist = np.unique(self.chaos_init_dist_data)

        for mass_ in in_mass:           #For every initial mass, we will plot a graph depending on rSMBH
            init_mass_idx = np.where((self.chaos_init_mass_data == mass_).all(1))[0] #Find indices where data files correspond to the correct initial masses
            init_mass_idx_stab = np.where((self.init_parti_m == mass_).all(1))[0]
            
            fig = plt.figure(figsize=(9, 14))
            ax1 = fig.add_subplot(211)
            ax2 = fig.add_subplot(212)

            tot_pop  = [ ]
            tot_popS = [ ]
            y_max = [ ]  
            iter = -1
            colourcycler = cycle(colour_picker())

            for dist_ in in_dist:
                iter += 1
                fin_parti = self.chaos_fin_parti_data[init_mass_idx] #Filter the needed data based on the data files satisfying condition
                stab_time = self.chaos_stab_time_data[init_mass_idx]
                init_dist = self.chaos_init_dist_data[init_mass_idx]

                ini_parti = self.ini_parti_data[init_mass_idx_stab]
                inj_event = self.inj_event_data[init_mass_idx_stab]
                merge_event = self.merge_no_data[init_mass_idx_stab]
                init_distS = self.initial_dist[init_mass_idx_stab]

                dist_arrays = np.where(init_dist == dist_) #Find indices corresponding to the correct distances. 
                fin_parti = fin_parti[dist_arrays]         #Here we split the data relative to their distances.
                stab_time = stab_time[dist_arrays]
                pop_size, pop_samples = np.unique(fin_parti, return_counts=True) #Count the number of unique final populations
                
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
                            surv_rate.append(100 * pop_samplesS[i]/(pop_samplesS[j]+pop_samples[i]))
                            break
                    if (no_surv):
                        surv_rate.append(0)
                        pop_samplesS.insert(j, 1)
                        pop_sizeS.insert(j, pop_size[i])

                for j in range(len(pop_sizeS)):
                    only_surv = True
                    for i in range(len(pop_size)):
                        if pop_sizeS[j] == pop_size[i]:
                            only_surv = False
                            break
                    if (only_surv):
                        surv_rate.insert(j, 100)
                pop_samplesS = np.array(pop_samplesS)
                pop_sizeS = np.array(pop_sizeS)

                colours = next(colourcycler)
                N_parti_avg = [ ]
                if len(pop_size) == 0:
                    pass
                else:
                    pop_id = np.argwhere(pop_size > 2)  #Only look at samples with N>2
                    pop_size = pop_size[pop_id]
                    pop_samples = pop_samples[pop_id]
                    tot_pop.append(max(pop_size))

                    pop_idS = np.argwhere(pop_sizeS > 2)
                    pop_sizeS = pop_sizeS[pop_idS]
                    pop_samplesS = pop_samplesS[pop_idS]

                    for pop_, samp_ in zip(pop_size, pop_samples):
                        N_parti = np.argwhere(fin_parti == pop_)         #For every N, we will gather their indices where they are in the pickle file
                        N_parti_avg.append(np.mean(stab_time[N_parti]))  #This enables us to compute their average between the whole data set                     
                    
                    y_max.append(max(N_parti_avg))
                    ax1.scatter(pop_size, N_parti_avg, color = colours, edgecolor = 'black', label = r'$r_{SMBH} =$'+str(dist_)+' pc')
                    ax2.scatter(tot_popS, surv_rate, color = colours, edgecolors = 'black')
                    for j, xpos in enumerate(pop_size):
                        if iter == 0:
                            ax1.text(xpos, -0.08*max(N_parti_avg), '# Ejec.\n'+'Set '+str(iter)+': '+str('{:.0f}'.format(pop_samples[j][0])), fontsize = 'xx-small', ha = 'center' )
                        else:
                            ax1.text(xpos, -0.08*max(N_parti_avg)*(1+0.6*iter), 'Set '+str(iter)+': '+str(pop_samples[j][0]), fontsize = 'xx-small', ha = 'center' )
                    for j, xpos in enumerate(tot_popS):
                        if iter == 0:
                            ax2.text(xpos, -0.08*max(surv_rate), '# Surv.\n'+'Set '+str(iter)+': '+str('{:.0f}'.format(pop_samplesS[j][0])), fontsize = 'xx-small', ha = 'center' )
                        else:
                            ax2.text(xpos, -0.08*max(surv_rate)*(1+0.6*iter), 'Set '+str(iter)+': '+str(pop_samplesS[j][0]), fontsize = 'xx-small', ha = 'center' )
                xints = [i for i in range(1+int(max(tot_pop)))]

            order = int('{:.0f}'.format(np.log10(mass_)))
            ax1.text(3.2, 1.125*max(y_max), r'$m_{i} \in$ '+str('{:.0f}'.format(mass_/10**order))+r'$\times 10^{}$'.format(order)+r'$M_\odot$')
            
            slope, intercept, r_value, p_value, std_err = stats.linregress(tot_popS, surv_rate)
            xtemp = np.linspace(2.6, max(tot_popS)+1)
            ytemp = [slope*xval_ + intercept for xval_ in xtemp]
            ax2.plot(xtemp, ytemp, color = 'black', zorder = 1)
            ax2.text(3.2, 10, r'Line of Best Fit: $S = $'+str('{:.2f}'.format(slope))+r'$N$')

            self.stab_plotter_logistics(ax1, ax2, tot_pop, xints, y_max)
            for ax_ in [ax1, ax2]:
                plot_ini.tickers(ax_)
            ax1.xaxis.labelpad = 25
            ax2.xaxis.labelpad = 25
            ax1.legend()
            plt.savefig('figures/chaotic_stab_time_equal_mass_'+str(mass_)+'_mean.pdf', dpi = 300, bbox_inches='tight')

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
                    ax.text(2.8, 0.96*(max(np.add(N_parti_avg, std))), r'$m_{IMBH}=$'+str('{:.0f}'.format(mass_/10**order))+r'$\times 10^{}$'.format(order)+r' $M_\odot$'+'\n'+r'$r_{SMBH}=$'+str(dist_)+' pc')
                         
                    self.error_plots(ax, xints, tot_pop, N_parti_avg, std)
                    plot_ini.tickers(ax)  
                    plt.savefig('figures/chaotic_stab_time_equal_mass_'+str(mass_)+'_err_dist_'+str(dist_)+'.pdf', dpi = 300, bbox_inches='tight')

    def steadytime_massdep_plotter(self):

        plot_ini = plotter_setup()

        in_dist = np.unique(self.chaos_init_dist_data)
        in_mass = np.unique(self.chaos_init_mass_data, axis=0)

        for dist_ in in_dist:
            init_dist_idx_chaos = np.where((self.chaos_init_dist_data == dist_))
            init_dist_idx_stab  = np.where((self.initial_dist == dist_))

            fig = plt.figure(figsize=(9, 14))
            ax1 = fig.add_subplot(211)
            ax2 = fig.add_subplot(212)

            tot_pop  = [ ]
            tot_popS = [ ]
            y_max = [ ]     
            iter = -1   
            coloursycler = cycle(colour_picker())

            for mass_ in in_mass:
                iter += 1

                init_mass = self.chaos_init_mass_data[init_dist_idx_chaos]
                init_mass = self.chaos_init_mass_data[init_dist_idx_chaos[0]]
                fin_parti = self.chaos_fin_parti_data[init_dist_idx_chaos[0]]
                stab_time = self.chaos_stab_time_data[init_dist_idx_chaos[0]]

                init_massS = self.init_parti_m[init_dist_idx_stab]
                ini_parti = self.ini_parti_data[init_dist_idx_stab[0]]
                inj_event = self.inj_event_data[init_dist_idx_stab[0]]
                merge_event = self.merge_no_data[init_dist_idx_stab[0]]

                mass_arrays = np.where((init_mass == mass_).all(1))[0]  #Indices with the correct mass column
                fin_parti = fin_parti[mass_arrays]
                stab_time = stab_time[mass_arrays]
                pop_size, pop_samples = np.unique(fin_parti, return_counts=True)

                mass_arrays = np.where((init_massS == mass_).all(1))[0]  #Indices with the correct mass column
                ini_parti = ini_parti[mass_arrays]
                inj_event = inj_event[mass_arrays]
                merge_event = merge_event[mass_arrays]
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
                            surv_rate.append(100 * pop_samplesS[i]/(pop_samplesS[j]+pop_samples[i]))
                            break
                    if (no_surv):
                        surv_rate.append(0)
                        pop_samplesS.insert(j, 1)
                        pop_sizeS.insert(j, pop_size[i])
                        
                for j in range(len(pop_sizeS)):
                    only_surv = True
                    for i in range(len(pop_size)):
                        if pop_sizeS[j] == pop_size[i]:
                            only_surv = False
                            break
                    if (only_surv):
                        surv_rate.insert(j, 100)
                pop_samplesS = np.array(pop_samplesS)
                pop_sizeS = np.array(pop_sizeS)

                colours = next(coloursycler)
                N_parti_avg = [ ]
                std = [ ]
                if len(pop_size) == 0:
                    pass
                else:
                    pop_id = np.argwhere(pop_size > 2)
                    pop_size = pop_size[pop_id]
                    pop_samples = pop_samples[pop_id]
                    tot_pop.append(max(pop_size))

                    pop_idS = np.argwhere(pop_sizeS > 2)
                    pop_sizeS = pop_sizeS[pop_idS]
                    pop_samplesS = pop_samplesS[pop_idS]

                    for pop_, samp_ in zip(pop_size, pop_samples):
                        N_parti = np.argwhere(fin_parti == pop_)
                        N_parti_avg.append(np.mean(stab_time[N_parti]))
                        std.append(np.std(stab_time[N_parti]))

                    y_max.append(max(N_parti_avg))
                    ax1.scatter(pop_size, N_parti_avg, color = colours, edgecolor = 'black',
                            label = r'$m_{i} \in$ ['+str(mass_[0])+', '+str(mass_[1])+r'] $M_\odot$')
                    ax2.scatter(tot_popS, surv_rate, color = colours, edgecolors = 'black', zorder = 2)

                    for j, xpos in enumerate(pop_size):
                        if iter == 0:
                            ax1.text(xpos, -0.08*max(N_parti_avg), '# Ejec.\n'+'Set '+str(iter)+': '+str('{:.0f}'.format(pop_samples[j][0])), fontsize = 'xx-small', ha = 'center' )
                        else:
                            ax1.text(xpos, -0.08*max(N_parti_avg)*(1+0.6*iter), 'Set '+str(iter)+': '+str(pop_samples[j][0]), fontsize = 'xx-small', ha = 'center' )
                    for j, xpos in enumerate(pop_sizeS):
                        if iter == 0:
                            ax2.text(xpos, -0.08*max(surv_rate), '# Surv.\n'+'Set '+str(iter)+': '+str('{:.0f}'.format(pop_samplesS[j][0])), fontsize = 'xx-small', ha = 'center' )
                        else:
                            ax2.text(xpos, -0.08*max(surv_rate)*(1+0.6*iter), 'Set '+str(iter)+': '+str(pop_samplesS[j][0]), fontsize = 'xx-small', ha = 'center' )
                xints = [i for i in range(1+int(max(tot_pop)))]

            slope, intercept, r_value, p_value, std_err = stats.linregress(tot_popS, surv_rate)
            xtemp = np.linspace(2.6, max(tot_popS)+1)
            ytemp = [slope*xval_ + intercept for xval_ in xtemp]
            ax2.plot(xtemp, ytemp, color = 'black', zorder = 1)
            ax2.text(3.2, 10, r'Line of Best Fit: $S = $'+str('{:.2f}'.format(slope))+r'$N$')

            ax1.xaxis.labelpad = 25
            ax2.xaxis.labelpad = 25
            ax1.text(2.8, 1.125*max(y_max), r'$r_{SMBH}=$'+str(dist_)+' pc')
            self.stab_plotter_logistics(ax1, ax2, tot_pop, xints, y_max)
            ax1.legend()
            for ax_ in [ax1, ax2]:
                plot_ini.tickers(ax_)
            plt.savefig('figures/chaotic_stab_time_equal_dist_'+str(dist_)+'_mean.pdf', dpi = 300, bbox_inches='tight')

        for dist_ in in_dist:
            for mass_ in in_mass: 
                y_max = []

                plt.clf()
                fig, ax = plt.subplots(figsize=(9,6))

                init_dist_idx = np.where((self.chaos_init_dist_data == dist_))
                init_mass = self.chaos_init_mass_data[init_dist_idx]
                fin_parti = self.chaos_fin_parti_data[init_dist_idx]
                stab_time = self.chaos_stab_time_data[init_dist_idx]

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

                    for j, xpos in enumerate(pop_size):
                        ax.text(xpos, -0.08*max(np.add(N_parti_avg, std)), '# Sim:\n'+str('{:.0f}'.format(pop_samples[j][0])), fontsize = 'xx-small', ha = 'center' )                   
                    ax.text(2.8, 0.96*(max(np.add(N_parti_avg, std))), r'$r_{SMBH}=$'+str(dist_)+' pc\n'+r'$m_{i} \in$ ['+str(mass_[0])+', '+str(mass_[1])+r'] $M_\odot$')
                    xints = [i for i in range(1+int(max(tot_pop)))]
                    
                    self.error_plots(ax, xints, tot_pop, N_parti_avg, std)
                    plot_ini.tickers(ax)
                    plt.savefig('figures/chaotic_stab_time_equal_dist_'+str(dist_)+'_err_mass_'+str(mass_)+'.pdf', dpi = 300, bbox_inches='tight')

gc_code = globular_cluster()
spatial_plotter(1.15*gc_code.gc_dist)
energy_plotter()
#animator(1.5 | units.parsec)

st = stability_time()
st.steadytime_massdep_plotter()
st.steadytime_distdep_plotter()
