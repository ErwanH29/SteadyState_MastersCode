from file_logistics import *
from spatial_plotters import *
from amuse.ext.galactic_potentials import MWpotentialBovy2015
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

class data_ext_files(object):
    """
    Class which extracts EJECTED data files
    """
    def __init__(self):

        self.IMBH_data, self.ejec_data = ejected_stat_extractor('data/Hermite/no_addition/chaotic_simulation/*', 'Hermite')
        self.IMBH_data_GRX, self.ejec_data_GRX = ejected_stat_extractor('data/GRX/no_addition/chaotic_simulation/*', 'GRX')

class KE_PE_plotters(object):
    """
    Class to plot KE vs. PE
    """

    def energy_data_extract(self, filter, file_IMBH, file_ejec):
        """
        Function to extract particle energy data
        
        Inputs:
        filter:     String 'B' or None which decides whether you are looking (B)elow a certain threshold or not
        file_IMBH:  File which holds the complete information of the simulation
        file_ejec:  File which holds the summarising information on the simulations
        """

        MW_code = MWpotentialBovy2015()

        ext = np.empty(len(file_IMBH))
        eyt = np.empty(len(file_IMBH))
        ezt = np.empty(len(file_IMBH))
        vesc = np.empty(len(file_IMBH))

        ejec_KEt = np.empty(len(file_IMBH))
        ejec_PEt = np.empty(len(file_IMBH))
        Nclose = np.empty(len(file_IMBH))
        ejected = np.empty(len(file_IMBH))

        print(file_IMBH)

        for i in range(len(file_IMBH)):
            ext[i], eyt[i], ezt[i], vesc[i], ejec_KEt[i], ejec_PEt[i], Nclose[i], \
            ejected[i] = ejected_extract_final(file_IMBH[i], file_ejec[i], 'E')
        ejec_PEt = np.asarray(ejec_PEt)
        ejec_KE = np.asarray(ejec_KEt[np.isfinite(ejec_KEt)])
        ejec_PE = np.asarray(ejec_PEt[np.isfinite(ejec_KEt)])
        ex = ext[np.isfinite(ext)]
        ey = eyt[np.isfinite(eyt)] 
        ez = ezt[np.isfinite(ezt)]

        PE_MW = ((1e3 | units.MSun) * MW_code.get_potential_at_point(0, 1 | units.pc, 1 | units.pc, 1 | units.pc)).value_in(units.J) / max(ejec_KE)

        if filter == 'B':
            ejected_KE = np.asarray(ejec_KE[(ejec_KE < 1e46) & (abs(ejec_PE) < 1e46)])
            ejected_PE = np.asarray(ejec_PE[(ejec_KE < 1e46) & (abs(ejec_PE) < 1e46)])
            final_pos = np.asarray([(i**2+j**2+z**2)**0.5 for i, j, z in zip(ex, ey, ez)])[(ejec_KE < 1e46) & (abs(ejec_PE) < 1e46)]
            ejected_KE /= max(ejected_KE)
            ejected_PE /= max(ejected_KE)

        else:
            ejected_KE = np.asarray([i/np.nanmax(ejec_KE) for i in ejec_KE])
            ejected_PE = np.asarray([i/np.nanmax(ejec_KE) for i in ejec_PE])
            final_pos = np.asarray([(i**2+j**2+z**2)**0.5 for i, j, z in zip(ex, ey, ez)])

        final_pos[np.isnan(final_pos)] = 0
        linex = np.linspace(-0.01*np.min(ejected_KE), 10)
        liney = -1*linex

        return ejected_KE, ejected_PE, final_pos, linex, liney, PE_MW

    def energy_plot_setup(self, KE, PE, linex, liney, MW_line, cdata, clabel, no_axis, save_file, hist):
        """
        Function to plot the required data
        
        Inputs:
        KE/PE:       The x, y data of the plot
        linex/liney: The line of best fit to distinguish bound from unbound
        cdata:       The data for which the colour code is based off
        no_axis:     The number of axis plotting (zoom or not)
        save_file:   String to ensure no overwritten plots
        hist:        String to state whether plotting histogram or other coloured variable
        """

        plot_ini = plotter_setup()
        plt.figure(figsize=(10, 4))
        if no_axis == 2:
            gs = gridspec.GridSpec(1, 2)
            ax1 = plt.subplot(gs[0,0])
            ax2 = plt.subplot(gs[0,1])
            ax1.set_ylabel(r'$E_P/E_{{K, {\rm{max}}}}$')
           
            for ax_ in [ax1, ax2]:
                ax_.plot(linex, liney, color = 'black', linestyle = '-.', zorder = 1)
                ax_.set_xlabel(r'$E_K/E_{{K, {\rm{max}}}}$')
                plot_ini.tickers(ax_, 'plot')
                
            if hist == 'Y':
                ax.set_title(r'Ejected Particle $K_E$ vs. $P_E$')

                bin2d_sim, xedges_s, yedges_s, image = ax2.hist2d(KE, PE, bins=(200, 200), range=([0,1.1],[1.1*min(PE),0]))
                bin2d_sim /= np.max(bin2d_sim)
                extent = [0, 1.1, 1.1*min(PE), 0]

                contours = ax2.imshow(np.log10(bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
                bin2d_sim, xedges_s, yedges_s, image = ax1.hist2d(KE, PE, bins=(80, 80), range=([0,0.35],[-0.35,0]))
                bin2d_sim /= np.max(bin2d_sim)
                extent = [0, 0.3, -0.3,0]

                contours = ax1.imshow(np.log10(bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
                cbar = plt.colorbar(contours, ax=ax2)
                cbar.set_label(label = r'$\log_{10}(N/N_{\rm{max}})$', labelpad=15)

                ax1.set_xlim(0,0.3)
                ax1.set_ylim(-0.3,0)
                ax2.set_xlim(0,max(KE))
                plt.savefig('figures/ejection_stats/KEPE_diagram_'+str(save_file)+'_.pdf', dpi=500, bbox_inches='tight')

            if hist == 'N':
                colour_axes = ax1.scatter(KE, PE, edgecolor = 'black', c=cdata, s = 20, zorder = 2)
                colour_axes = ax2.scatter(KE, PE, edgecolor = 'black', c=cdata, s = 20, zorder = 3)
                plt.colorbar(colour_axes, ax=ax2, label = clabel)
                ax2.set_xlim(0,1.05)
                ax2.set_ylim(-0.1, 0)
                ax1.set_xlim(0,0.1)
                ax1.set_ylim(-0.1,0)
                
            return
        
        else:
            if hist == 'N':
                fig, ax = plt.subplots()
                ax.set_title(r'Ejected Particle $K_E$ vs. $P_E$')
                ax.set_xlabel(r'$E_K/E_{{K, \rm{max}}}$')
                ax.set_ylabel(r'$E_{{P}}/E_{{K, \rm{max}}}$')

                MW_linex = np.linspace(abs(MW_line), 1.05)
                MW_liney = [MW_line for i in range(len(MW_linex))]

                ax.plot(linex, liney, color = 'black', linestyle = '-.', zorder = 1)
                colour_axes = ax.scatter(KE, PE, c = cdata, edgecolors = 'black', s = 20, zorder = 3)
                ax.plot(MW_linex, MW_liney, linestyle = '--', color = 'black', zorder = 2)
                plt.colorbar(colour_axes, ax = ax, label = r'$r_{\rm{GC}}$ [pc]')

                ax.set_xlim(0,1.05)
                ax.set_ylim(-1.05,0)
                plot_ini.tickers(ax, 'hist')
                plt.savefig('figures/ejection_stats/KEPE_histogram_'+str(save_file)+'.pdf', dpi=300, bbox_inches='tight')
                plt.clf()
                
            else: 
                fig, ax = plt.subplots()

                bin2d_sim, xedges_s, yedges_s, image = ax.hist2d(KE, PE, bins=(100,300), range=([np.min(KE), 2],[1.1*min(PE),0.1]))
                bin2d_sim /= np.max(bin2d_sim)
                extent = [0, 1.1, 1.1*min(PE),0]

                ax.set_title(r'Ejected Particle \n $K_E$ vs. $P_E$')
                ax.set_xlim(1e-6, 1.1)
                ax.plot(linex, liney, color = 'black', linestyle = '-.', zorder = 1)
                ax.set_ylim(-1, -10**-6)
                ax.set_xlabel(r'$E_K/E_{{K, {}}}$'.format(str('max')))
                ax.set_ylabel(r'$E_P/E_{{K, {}}}$'.format(str('max')))
                contours = ax.imshow(np.log10(bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
                ax.plot(linex, liney, color = 'black', linestyle = '-.')
                cbar = plt.colorbar(contours, ax=ax)
                cbar.set_label(r"$\log_{10}(N/N_{max})$", labelpad=15)
                
                ax.set_yscale('symlog')
                plot_ini.tickers(ax, 'plot')
                plt.savefig('figures/ejection_stats/KEPE_diagram_'+str(save_file)+'.pdf', dpi=300, bbox_inches='tight')
                plt.clf()
            return

    def KEPE_plotter(self):
            """
            Function to plot the KE vs. PE plot coloured with final positions or with histogram
            """

            data = data_ext_files()

            save_file = ['fpos', 'fpos_crop']
            data_filt = [None, 'B']

            for i in range(1):
                axis = i+1
                ejected_KE, ejected_PE, final_pos, unbounded_x, unbounded_y, line_MW = self.energy_data_extract(data_filt[i], data.IMBH_data, data.ejec_data)
                self.energy_plot_setup(ejected_KE, ejected_PE, unbounded_x, unbounded_y, line_MW, final_pos, 'Final Distance to Core [pc]', axis, save_file[i], 'N')

class vejection(object):
    """
    Class to plot the velocity ejection statistics
    """

    def vejec_plotter(self):
        """
        Function to plot how the total population/ mass of the system influences the ejection velocity
        """
        plot_ini = plotter_setup()
        data = data_ext_files()
        plot_ini = plotter_setup()
        MW_code = MWpotentialBovy2015()

        ex = [[ ], [ ]]
        ey = [[ ], [ ]]
        ez = [[ ], [ ]]
        vesc = [[ ], [ ]]

        ejec_KE = [[ ], [ ]]
        ejec_PE = [[ ], [ ]]
        Nclose = [[ ], [ ]]

        tot_pop = [[ ], [ ]]
        ejected = [[ ], [ ]]
        surv_time = [[ ], [ ]]
        avg_surv = [[ ], [ ]]
        avg_vel = [[ ], [ ]]

        samples = [[ ], [ ]]
        pops = [[ ], [ ]]
        IMBH_data = [data.IMBH_data, data.IMBH_data_GRX]
        ejec_data = [data.ejec_data, data.ejec_data_GRX]

        for int_ in range(1):
            for i in range(len(IMBH_data[int_])):
                ex_val, ey_val, ez_val, vesc_val, ejec_KE_val, ejec_PE_val, Nclose_val, \
                ejected_val = ejected_extract_final(IMBH_data[int_][i], ejec_data[int_][i], 'E')

                ex[int_].append(ex_val)
                ey[int_].append(ey_val)
                ez[int_].append(ez_val)
                vesc[int_].append(vesc_val)
                ejec_KE[int_].append(ejec_KE_val)
                ejec_PE[int_].append(ejec_PE_val)
                Nclose[int_].append(Nclose_val)
                ejected[int_].append(ejected_val)

                vals_df = ejec_data[int_][i].iloc[0]
                tot_pop[int_].append(vals_df[6])
                surv_time[int_].append(vals_df[-2].value_in(units.Myr))
            vesc[int_] = np.asarray([i for i in vesc[int_]])

            in_pop = np.unique(tot_pop[int_])
            avg_vel_t = np.empty(len(in_pop))
            avg_surv_t = np.empty(len(in_pop))

            iter = -1
            for pop_ in in_pop:
                iter += 1
                pops[int_].append(pop_)
                indices = np.where((tot_pop[int_] == pop_))[0]
                samples[int_].append(len(vesc[int_][indices]))
                avg_vel_t[iter] = np.mean([vesc[int_][i] for i in indices])
                avg_surv_t[iter] = np.mean([surv_time[int_][i] for i in indices])
            avg_vel[int_] = [i for i in avg_vel_t]
            avg_surv[int_] = [i for i in avg_surv_t]
        vesc_MW = (np.sqrt(2)*MW_code.circular_velocity(0.1 | units.parsec) + np.sqrt(2*constants.G*(4e6 | units.MSun)/(0.1 | units.parsec))).value_in(units.kms)
        
        fig = plt.figure(figsize=(15, 6))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax_ = [ax1, ax2]
        ax_title = ['Hermite', 'GRX']    
        ymin = min(avg_vel[0])
        ymax = max(avg_vel[0])
        """ymax = max(max(avg_vel[0]), max(avg_vel[1]))
        norm_min = np.log10(min(min(avg_surv[0]), min(avg_surv[1])))
        norm_max = np.log10(max(max(avg_surv[0]), max(avg_surv[1])))"""
        norm_min = np.log10(min(avg_surv[0]))#, min(avg_surv[1])))
        norm_max = np.log10(max(avg_surv[0]))#, max(avg_surv[1])))
        normalise = plt.Normalize(norm_min, norm_max)
        for int_ in range(1):
            colour_axes = ax_[int_].scatter(pops[int_], avg_vel[int_], edgecolors='black', c = np.log10(avg_surv[int_]), norm = normalise, zorder = 3)
            ax_[int_].set_ylabel(r'$\langle v_{ej} \rangle$ [km/s]')
            ax_[int_].axhline(vesc_MW, color = 'black', linestyle = ':')
            ax_[int_].text(15, 660, r'$v_{\rm{esc, MW}}$')
            ax_[int_].set_ylim(0.9*ymin, 1.1*ymax)
            ax_[int_].set_title(ax_title[int_])
            ax_[int_].xaxis.labelpad = 30
            for j, xpos in enumerate(pops[int_]):
                ax_[int_].text(pops[int_][j], 0.6*ymin, str(ax_title[int_])+'\n'+str('{:.0f}'.format(samples[int_][j])), fontsize = 'xx-small', ha = 'center' )
            plot_ini.tickers_pop(ax_[int_], pops[int_])
        plt.colorbar(colour_axes, ax=ax2, label = r'$\log_{10}\langle t_{ej}\rangle$ [Myr]')
        plt.savefig('figures/ejection_stats/mean_vej.pdf', dpi=300, bbox_inches='tight')

        fig = plt.figure(figsize=(15, 6))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax_ = [ax1, ax2]
        ax_title = ['Ejection Velocity Histogram \n Hermite', '\n GRX']
        colours = ['red', 'blue']
        #xmax = max(max(vesc[0]), max(vesc[1]))
        xmax = (max(vesc[0]))
        
        for int_ in range(1):
            n, bins, patches = ax_[int_].hist(vesc[int_], 20, histtype = 'step', color=colours[int_])
            n, bins, patches = ax_[int_].hist(vesc[int_], 20, color=colours[int_], alpha = 0.4)
            
            ax_[int_].axvline(vesc_MW, linestyle = ':', color = 'black')
            ax_[int_].text(vesc_MW*(1+0.05), 0.9*max(n), r'$v_{esc, MW}$', rotation = 270, horizontalalignment = 'center')
            ax_[int_].set_xlabel(r'$v_{ejec}$ [km s$^{-1}$]')
            ax_[int_].set_ylabel(r'Frequency')
            ax_[int_].set_xlim(0,1.01*xmax)
            plot_ini.tickers(ax_[int_], 'plot')
        plt.savefig('figures/ejection_stats/vejection_histogram.pdf', dpi = 300, bbox_inches='tight')

class event_tracker(object):
    """
    Class to take stats of ejection vs. mergers
    """
    def __init__(self):

        plot_ini = plotter_setup()
        chaos_data = glob.glob('data/Hermite/no_addition/chaotic_simulation/*')
        chaos_data_GRX = glob.glob('data/GRX/no_addition/chaotic_simulation/*')
        chaos_data = [natsort.natsorted(chaos_data), natsort.natsorted(chaos_data_GRX)]

        init_pop = [[ ], [ ]]
        merger = [[ ], [ ]]
        in_pop = [[ ], [ ]]
        frac_merge = [[ ], [ ]]
        colours = ['red', 'blue']
        labels = ['Hermite', 'GRX']

        fig, ax = plt.subplots()
        ax.set_ylabel(r'$N_{\rm{merge}}/N_{\rm{sim}}$')
        ax.set_ylim(0,1.05)
        for int_ in range(2):
            for file_ in range(len(chaos_data[int_])):
                with open(chaos_data[int_][file_], 'rb') as input_file:
                    data = pkl.load(input_file)
                    init_pop[int_].append(len(data.iloc[0][8])+1)
                    merger[int_].append(data.iloc[0][10])

            in_pop[int_] = np.unique(init_pop[int_])
            iter = -1
            for pop_ in in_pop[int_]:
                iter +=1
                indices = np.where((init_pop[int_] == pop_))[0]
                temp_frac = [merger[int_][i] for i in indices]
                frac_merge[int_].append(np.mean(temp_frac))

            ax.scatter(in_pop[int_], frac_merge[int_], color = colours[int_], edgecolors = 'black', label = labels[int_])
        plot_ini.tickers_pop(ax, in_pop[0])
        plt.legend()
        plt.savefig('figures/ejection_stats/SMBH_merge_fraction.pdf', dpi=300, bbox_inches='tight')


print('...ejection_stat_plotters...')
cst = event_tracker()
cst = vejection()
cst.vejec_plotter()
cst = KE_PE_plotters()
cst.KEPE_plotter()