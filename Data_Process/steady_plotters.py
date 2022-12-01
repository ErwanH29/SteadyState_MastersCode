from amuse.lab import *
from file_logistics import *
from spatial_plotters import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
from scipy.optimize import OptimizeWarning

class stability_plotters(object):
    """
    Class to extract data and setup plots for simulations where IMBH particles are added over time
    """
        
    def chaos_extract(self, dirC):
        """
        Function to extract data from simulations who end in ejection
        """

        warnings.filterwarnings("ignore", category=OptimizeWarning) 
        self.chaos_ini_parti_data, self.chaos_fparti_data, self.chaos_number_mergers, self.chaos_cumulative_mm, self.chaos_simulated_end, \
        self.chaos_ejected_parti, self.chaos_stab_time_data, self.chaos_init_dist_data, self.chaos_init_mass_data, \
        self.chaos_inj_mass_data, self.chaos_eje_mass_data = stats_chaos_extractor(dirC)  

        return self.chaos_ini_parti_data, self.chaos_fparti_data, self.chaos_number_mergers, self.chaos_cumulative_mm, \
               self.chaos_simulated_end, self.chaos_ejected_parti, self.chaos_stab_time_data, self.chaos_init_dist_data, \
               self.chaos_init_mass_data, self.chaos_inj_mass_data, self.chaos_eje_mass_data

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
        pop, psamples = np.unique(filt_finparti, return_counts = True)

        return pop, psamples, filt_finparti, filt_stabtime

    def stable_extract(self, dirS):
        """
        Function to extract data  from simulations who end with stable state
        """
        self.ini_parti_data, self.inj_event_data, self.merge_no_data, self.mergerm_data, self.simulated_end, \
        self.initial_dist, self.cluster_rad, self.init_parti_m, self.trelax_data = stats_stable_extractor(dirS)

        return self.ini_parti_data, self.inj_event_data, self.merge_no_data, self.mergerm_data, self.simulated_end, \
               self.initial_dist, self.cluster_rad, self.init_parti_m, self.trelax_data
    

    def spread_steady_plotter(self, in_dist, in_mass, idist_data, imass_data, fparti_data, stab_time_data, integrator):
        """
        Function to plot the spread in the stability time
        
        Inputs:
        in_dist:         The initial distance of the cluster from the central SMBH
        in_mass:         The mass of the IMBH
        idist_data:      The initial distance of the cluster (array)
        imass_data:      The mass of the IMBH (array)
        fparti_data:     Data of the particles present at end time
        stab_time_data:  Array consisting of total simulation time
        integrator:      Integrator (Hermite or GRX) used
        """
        plot_ini = plotter_setup()

        for dist_ in in_dist:
            for mass_ in in_mass: 
                y_max = []

                plt.clf()
                fig, ax = plt.subplots(figsize=(8,6))

                init_dist_idx = np.where((idist_data == dist_))
                init_mass = imass_data[init_dist_idx]
                fparti = fparti_data[init_dist_idx]
                stab_time = stab_time_data[init_dist_idx]

                tot_pop = [ ]
                N_parti_avg = [ ]
                std = [ ]
                mass_arrays = np.where((init_mass == mass_).all(1))[0]  #Indices with the correct mass column

                fparti = fparti[mass_arrays]
                stab_time = stab_time[mass_arrays]
                pop, psamples = np.unique(fparti, return_counts=True)
                pop_id = np.argwhere(pop > 2)
                pop = pop[pop_id]

                if len(pop) == 0:
                    pass
                else:
                    psamples = psamples[pop_id]
                    tot_pop.append(max(pop))

                    for pop_, samp_ in zip(pop, psamples):
                        N_parti = np.argwhere(fparti == pop_)
                        N_parti_avg.append(np.mean(stab_time[N_parti]))
                        std.append(np.std(stab_time[N_parti]))
                    y_max.append(max(np.add(N_parti_avg, std)))

                    ax.errorbar(pop, (N_parti_avg), color = 'black', yerr=std, fmt = 'o')
                    ax.scatter(pop, (np.add(N_parti_avg, std)), marker = '_', color = 'black')
                    ax.scatter(pop, (np.subtract(N_parti_avg, std)), marker = '_', color = 'black')
                    ax.scatter(pop, (N_parti_avg), color = 'black')
                    
                    ax.text(82, (87), r'$r_{SMBH}=$'+str(dist_)+' pc\n'+r'$m_{i} =$ '+str(mass_[0])+r' $M_\odot$')
                    ax.set_xlim(5, 105)
                    ax.set_ylim(10**-3, 120)
                    ax.set_ylabel(r'$\log_{10} t_{\rm{eject}}$ [Myr]')
                    ax.set_title(r'Spread in Stability Time')
                    plot_ini.tickers_pop(ax, pop)
                    plt.savefig('figures/steady_time/spread_steady_time'+str(integrator)+'_dist_'+str(dist_)+'.pdf', dpi = 300, bbox_inches='tight')

    def overall_steady_plotter(self):
        """
        Function to plot stability time for constant distances
        """

        def log_slope(slope):
            rclust = 0.1 | units.pc
            Mclust = 4.1*10**6 | units.MSun
            slope = slope * (1 | units.Myr)
            trlx_coeff = (rclust/(constants.G*Mclust)**1/3)
            return np.log(abs(slope.value_in(units.s)))/np.log(trlx_coeff.value_in(units.s**2/units.m**2))
            
        def log_fit(xval, slope, alpha, yint):
            return (slope)/(np.log(alpha*xval)) + yint

        plot_ini = plotter_setup()

        dirH = 'data/Hermite/no_addition/chaotic_simulation/*'
        dirG = 'data/GRX/no_addition/chaotic_simulation/*'

        chaos_ini_parti_data, chaos_fparti_data, chaos_number_mergers, chaos_cumulative_mm, chaos_simulated_end, \
        chaos_ejected_parti, chaos_stab_time_data, chaos_init_dist_data, chaos_init_mass_data, \
        chaos_inj_mass_data, chaos_eje_mass_data = self.chaos_extract(dirH)

        chaos_ini_parti_dataG, chaos_fparti_dataG, chaos_number_mergersG, chaos_cumulative_mmG, chaos_simulated_endG, \
        chaos_ejected_partiG, chaos_stab_time_dataG, chaos_init_dist_dataG, chaos_init_mass_dataG, \
        chaos_inj_mass_dataG, chaos_eje_mass_dataG = self.chaos_extract(dirG)

        in_dist = np.unique(chaos_init_dist_data)
        in_mass = np.unique(chaos_init_mass_data, axis=0)

        for dist_ in in_dist:
            init_dist_idx_chaos = np.where((np.asarray(chaos_init_dist_data) == dist_))
            init_dist_idx_chaosG = np.where((np.asarray(chaos_init_dist_dataG) == dist_))

            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111)

            tot_pop  = [ ]   
            iter = -1   
            for mass_ in in_mass:
                iter += 1
                pop, psamples, fparti, stab_time = self.index_extractor(chaos_init_mass_data, chaos_fparti_data, 
                                                                        chaos_stab_time_data, init_dist_idx_chaos, mass_)
                popG, psamplesG, fpartiG, stab_timeG = self.index_extractor(chaos_init_mass_dataG, chaos_fparti_dataG, 
                                                                            chaos_stab_time_dataG, init_dist_idx_chaosG, mass_)

                N_parti_avg = [ ]
                N_parti_std = [ ]
                pop_id = np.argwhere(pop > 2)
                pop = pop[pop_id]
                psamples = psamples[pop_id]

                N_parti_avgG = [ ]
                N_parti_stdG = [ ]
                pop_idG = np.argwhere(popG > 2)
                popG = popG[pop_idG]
                psamplesG = psamplesG[pop_idG]

                tot_pop.append(max(pop))

                full_simul_H = []
                for pop_, samp_ in zip(pop, psamples):
                    N_parti = np.argwhere(fparti == pop_)
                    N_parti_avg.append(np.mean(stab_time[N_parti]))
                    N_parti_std.append(np.std(stab_time[N_parti]))
                    idx = np.where(stab_time[N_parti] == 100)[0]
                    ratio = len(idx)/len(stab_time[N_parti])
                    full_simul_H.append(ratio)
                N_parti_avg = np.asarray(N_parti_avg)

                full_simul_G = []
                for pop_, samp_ in zip(popG, psamplesG):
                    N_parti = np.argwhere(fpartiG == pop_)
                    N_parti_avgG.append(np.mean(stab_timeG[N_parti]))
                    N_parti_stdG.append(np.std(stab_timeG[N_parti]))
                    idx = np.where(stab_timeG[N_parti] == 100)[0]
                    ratio = len(idx)/len(stab_timeG[N_parti])
                    full_simul_G.append(ratio)
                N_parti_avgG = np.asarray(N_parti_avgG)

                """
                ax.scatter(pop, np.log10(N_parti_avg), color = 'red', edgecolor = 'black', zorder = 3,
                           label = r'Hermite')                
                ax.scatter(popG, np.log10(N_parti_avgG), edgecolor = 'black', color = 'blue', 
                           zorder = 3, label = r'Hermite GRX')  """  
                           
                ax.errorbar(pop, N_parti_avg, color = 'red', yerr=N_parti_std, fmt = 'o', label = 'Hermite')
                ax.scatter(pop, N_parti_avg, color = 'black', s= 60)
                ax.scatter(pop, N_parti_avg+N_parti_std, marker = '_', color = 'red')
                ax.scatter(pop, N_parti_avg-N_parti_std, marker = '_', color = 'red')
                ax.errorbar(popG, N_parti_avgG, color = 'blue', yerr=N_parti_stdG, fmt = 'o', label = 'GRX')
                ax.scatter(popG, N_parti_avgG, color = 'black', s= 60)
                ax.scatter(popG, N_parti_avgG+N_parti_stdG, marker = '_', color = 'blue')
                ax.scatter(popG, N_parti_avgG-N_parti_stdG, marker = '_', color = 'blue')

                ax.set_ylabel(r'$\log_{10} t_{\rm{surv}}$ [Myr]')   

            for j, xpos in enumerate(pop):
                ax.text(pop[j][0], 3.25*10**-1, '# Ejec.\n'+'Hermite: '+str('{:.0f}'.format(psamples[j][0])), fontsize = 'xx-small', ha = 'center' )
            for j, xpos in enumerate(popG):
                ax.text(popG[j][0], 2.9*10**-1, 'GRX: '+str('{:.0f}'.format(psamplesG[j][0])), fontsize = 'xx-small', ha = 'center' )

            pop = np.array([float(i) for i in pop])
            N_parti_avg = np.array([float(i) for i in N_parti_avg])   
            plot_ini.tickers_pop(ax, pop)
            ax.set_xlim(5,105)


            ymax = [N_parti_avg[0]+N_parti_std[0], 14.558, 0]
            ymin = [N_parti_avg[0]-N_parti_std[0], 10.94627273, 14.558, N_parti_avg[-1]+N_parti_std[-1]]
            p0 = (1.5,  0.9, 10)
            params, cv = scipy.optimize.curve_fit(log_fit, pop, (N_parti_avg), p0, maxfev = 2000)
            params_max, cv = scipy.optimize.curve_fit(log_fit, [10, 50, 100], ymax, p0, maxfev = 2000)
            params_min, cv = scipy.optimize.curve_fit(log_fit, [10, 40, 50, 100], ymin, p0, maxfev = 2000)

            slope, alpha, intercept = params
            slope_max, alpha_max, intercept_max = params_max
            slope_min, alpha_min, intercept_min = params_min
            
            """paramsG, cv = scipy.optimize.curve_fit(log_fit, pop, N_parti_avg, p0) #TO CHANGE WHEN GRX
            slopeG, interceptG = paramsG"""

            x_arr = np.linspace(10,100)
            y_arr =([10**-1.64*np.log(alpha*i)**-1 for i in x_arr])

            slope_str = str('{:.2f}'.format(slope))
            alpha_str = str('{:.2f}'.format(alpha))
            xtemp = np.linspace(5, 105, 1000)
            y_bf = [log_fit(i, slope, alpha, intercept) for i in xtemp]
            y_bf1 = [log_fit(i, slope_min, alpha_min, intercept_min) for i in xtemp]
            y_bf2 = [log_fit(i, slope_max, alpha_max, intercept_max) for i in xtemp]
            ax.plot(xtemp, (y_bf), zorder = 1, color = 'black', ls = '-.')
            #ax.plot(xtemp, (y_bf1), zorder = 1, color = 'black', ls = ':')
            #ax.plot(xtemp, (y_bf2), zorder = 1, color = 'black', ls = ':')
            ax.plot(x_arr, (y_arr))
            ax.text(77, 52, r'$t_{{\rm surv}} \approx (\frac{{{}}}{{\ln{{{}}}N}}$'.format(slope_str, alpha_str)+') Myr')
            ax.set_ylim(5*10**-1, 1.25*(max(N_parti_avg+N_parti_std)))
            ax.set_yscale('log')
            ax.legend()
            ax.xaxis.labelpad = 30
            plt.savefig('figures/steady_time/const_population_stab_time_equal_dist_'+str(dist_)+'_mean.pdf', dpi = 300, bbox_inches='tight')

            errs_upper = log_slope(slope_max) - log_slope(slope)
            erra_upper = alpha_max - alpha
            errs_lower = log_slope(slope_min) - log_slope(slope)
            erra_lower = alpha_min - alpha
            with open('figures/steady_time/Sim_summary.txt', 'w') as file:
                file.write('For Hermite, # of full simulations per population: '+str(pop.flatten())+str(full_simul_H))
                file.write('\nThe slope of the full curve goes as:               '+str(log_slope(slope)))
                file.write('\nwith errors:                                       +'+str(errs_upper)+' '+str(errs_lower))
                file.write('\nThe logarithmic prefactor goes as:                 '+str(alpha))
                file.write('\nwith errors:                                       '+str(erra_upper)+' '+str(erra_lower))
                file.write('\n\nFor GRX, # of full simulations per population: '+str(pop.flatten())+str(full_simul_G))
                #file.write('\nand the slope of the curve goes as:                '+str(np.log(slopeG)*np.exp(1)/10))

        self.spread_steady_plotter(in_dist, in_mass, chaos_init_dist_data,
                                   chaos_init_mass_data, chaos_fparti_data,
                                   chaos_stab_time_data, 'Hermite')
        self.spread_steady_plotter(in_dist, in_mass, chaos_init_dist_dataG,
                                   chaos_init_mass_dataG, chaos_fparti_dataG,
                                   chaos_stab_time_dataG, 'GRX')

print('...steady_plotter...')
cst = stability_plotters()
cst.overall_steady_plotter()