from amuse.ext.galactic_potentials import MWpotentialBovy2015
from file_logistics import *
from spatial_plotters import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

class data_ext_files(object):
    """
    Class which extracts EJECTED data files
    """
    def __init__(self):

        self.IMBH_data, self.ejec_data = ejected_stat_extractor('data/Hermite/no_addition/chaotic_simulation/*')
        #self.IMBH_data_GRX, ejec_data_GRX = ejected_stat_extractor('data/Hermite/no_addition/chaotic_simulation/*')

        self.ex = np.empty((len(self.IMBH_data)))
        self.ey = np.empty((len(self.IMBH_data)))
        self.ez = np.empty((len(self.IMBH_data)))
        self.vesc = np.empty((len(self.IMBH_data)))

        self.ejec_KE = np.empty((len(self.IMBH_data)))
        self.ejec_PE = np.empty((len(self.IMBH_data)))
        self.Nclose = np.empty((len(self.IMBH_data)))

        self.tot_mass = np.empty((len(self.IMBH_data)))
        self.tot_pop = np.empty((len(self.IMBH_data)))
        self.surv_time = np.empty((len(self.IMBH_data)))
        self.ejected = np.empty(len(self.IMBH_data))

        """self.ex_GRX = np.empty((len(self.IMBH_data)))
        self.ey_GRX = np.empty((len(self.IMBH_data)))
        self.ez_GRX = np.empty((len(self.IMBH_data)))
        self.vesc_GRX = np.empty((len(self.IMBH_data)))

        self.ejec_KE_GRX = np.empty((len(self.IMBH_data)))
        self.ejec_PE_GRX = np.empty((len(self.IMBH_data)))
        self.Nclose_GRX = np.empty((len(self.IMBH_data)))

        self.tot_mass_GRX = np.empty((len(self.IMBH_data)))
        self.tot_pop_GRX = np.empty((len(self.IMBH_data)))
        self.surv_time_GRX = np.empty((len(self.IMBH_data)))
        self.ejected_GRX = np.empty(len(self.IMBH_data))"""

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

        ext = np.empty((len(file_IMBH)))
        eyt = np.empty((len(file_IMBH)))
        ezt = np.empty((len(file_IMBH)))
        vesc = np.empty((len(file_IMBH)))

        ejec_KEt = np.empty((len(file_IMBH)))
        ejec_PEt = np.empty((len(file_IMBH)))
        Nclose = np.empty((len(file_IMBH)))
        ejected = np.empty(len(file_IMBH))

        for i in range(len(file_IMBH)):
            ext[i], eyt[i], ezt[i], vesc[i], ejec_KEt[i], ejec_PEt[i], Nclose[i], \
            ejected[i] = ejected_extract_final(file_IMBH[i], file_ejec[i], 'E') #Change S to E
        ejec_PEt = np.asarray(ejec_PEt)
        #ejec_PE = np.asarray([PE_ + ((1000*1|units.MSun) * MW_code.get_potential_at_point(0, x * 1 | units.pc, y * 1 | units.pc, z * 1 | units.pc)).value_in(units.J) for PE_, x, y, z in zip(ejec_PE, ex, ey, ez)])
        ejec_KE = ejec_KEt[np.isfinite(ejec_KEt)]
        ejec_PE = ejec_PEt[np.isfinite(ejec_KEt)]
        ex = ext[np.isfinite(ext)] ; ey = eyt[np.isfinite(eyt)] ; ez = ezt[np.isfinite(ezt)]
        PE_MW = ((1e3 | units.MSun) * MW_code.get_potential_at_point(0, 1 | units.pc, 1 | units.pc, 1 | units.pc)).value_in(units.J) / max(ejec_KE)

        if filter == 'B':
            idx = np.where(ejec_KE < 1e46)                   # Change to 1e45 when more data
            ejected_KE1 = np.asarray([i for i in ejec_KE[idx]])
            ejected_PE1 = np.asarray([i for i in ejec_PE[idx]])
            final_pos = np.asarray([(i**2+j**2+z**2)**0.5 for i, j, z in zip(ex, ey, ez)])[idx]

            idx = np.where(abs(ejected_PE1) < 1e46)
            ejected_KE = np.asarray([i/max(ejected_KE1[idx]) for i in ejected_KE1[idx]])
            ejected_PE = np.asarray([i/max(ejected_KE1[idx]) for i in ejected_PE1[idx]])
            final_pos = np.asarray(final_pos[idx])

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
            ax1.set_ylabel(r'$E_P/E_{{K, {}}}$'.format(str('max')))
           
            for ax_ in [ax1, ax2]:
                ax_.plot(linex, liney, color = 'black', linestyle = '-.', zorder = 1)
                ax_.set_xlabel(r'$E_K/E_{{K, {}}}$'.format(str('max')))
                plot_ini.tickers(ax_, 'plot')
                
            if hist == 'Y':
                bin2d_sim, xedges_s, yedges_s, image = ax2.hist2d(KE, PE, bins=(200, 200), range=([0,1.1],[1.1*min(PE),0]))
                bin2d_sim /= np.max(bin2d_sim)
                extent = [0, 1.1, 1.1*min(PE), 0]

                ax.set_title(r'Ejected Particle $K_E$ vs. $P_E$')
                contours = ax2.imshow(np.log10(bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
                #ax2.scatter(KE, PE, s = 15, color = 'red', edgecolors = 'black')
                #ax2.set_ylim(m(PE),0)
                ax2.set_xlim(0,max(KE))
                bin2d_sim, xedges_s, yedges_s, image = ax1.hist2d(KE, PE, bins=(80, 80), range=([0,0.35],[-0.35,0]))
                bin2d_sim /= np.max(bin2d_sim)
                extent = [0, 0.3, -0.3,0]
                contours = ax1.imshow(np.log10(bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
                #ax1.scatter(KE, PE, s = 15, color = 'red', edgecolors = 'black')
                cbar = plt.colorbar(contours, ax=ax2)
                cbar.set_label(label = r'$\log_{10}(N/N_{max})$',rotation=270,labelpad=15)

                ax1.set_xlim(0,0.3)
                ax1.set_ylim(-0.3,0)
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

                MW_linex = np.linspace(abs(MW_line), 1.05)
                MW_liney = [MW_line for i in range(len(MW_linex))]

                fig, ax = plt.subplots()
                ax.set_title(r'Ejected Particle $K_E$ vs. $P_E$')
                ax.plot(linex, liney, color = 'black', linestyle = '-.', zorder = 1)
                colour_axes = ax.scatter(KE, PE, c = cdata, edgecolors = 'black', s = 20, zorder = 3)
                ax.set_xlim(0,1.05)
                ax.set_ylim(-1.05,0)
                ax.plot(MW_linex, MW_liney, linestyle = '--', color = 'black', zorder = 2)
                plt.colorbar(colour_axes, ax = ax, label = r'Final Distance to Core [pc]')
                ax.set_xlabel(r'$E_K/E_{{K, {}}}$'.format(str('max')))
                ax.set_ylabel(r'$E_{{P, MW}}/E_{{K, {}}}$'.format(str('max')))
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
                #ax.set_yscale('symlog')
                #ax.set_xscale('log')
                contours = ax.imshow(np.log10(bin2d_sim), extent = extent, aspect='auto', origin = 'upper')
                ax.plot(linex, liney, color = 'black', linestyle = '-.')
                #ax.scatter(KE, PE, s=15, color = 'red', edgecolors ='black')
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
            save_file_hist = ['histogram', 'histogram_crop']
            data_filt = [None, 'B']

            for i in range(2):
                axis = i+1
                ejected_KE, ejected_PE, final_pos, unbounded_x, unbounded_y, line_MW = self.energy_data_extract(data_filt[i], data.IMBH_data, data.ejec_data)
                self.energy_plot_setup(ejected_KE, ejected_PE, unbounded_x, unbounded_y, line_MW, final_pos, 'Final Distance to Core [pc]', axis, save_file[i], 'N')  
                #self.finalpos_energy_plotter(ejected_KE, ejected_PE, unbounded_x, unbounded_y, final_pos, r'$\log_{10}n$', i, save_file_hist[i], 'Y')  

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

        for i in range(len(data.ejec_data)):
            data.ex[i], data.ey[i], data.ez[i], data.vesc[i], data.ejec_KE[i], data.ejec_PE[i], data.Nclose[i], \
            data.ejected[i] = ejected_extract_final(data.IMBH_data[i], data.ejec_data[i], 'E')
            vals_df = data.ejec_data[i].iloc[0]
            data.tot_pop[i] = vals_df[6]
            data.surv_time[i] = vals_df[-2].value_in(units.Myr)

        in_pop = np.unique(data.tot_pop)
        avg_vel = np.empty(len(in_pop))
        avg_surv = np.empty((len(in_pop)))

        iter = -1
        for pop_ in in_pop:
            iter += 1
            indices = np.where((data.tot_pop == pop_))[0]
            temp_escv = [data.vesc[i] for i in indices]
            temp_surv = [data.surv_time[i] for i in indices]

            avg_vel[iter] = np.mean(temp_escv)
            avg_surv[iter] = np.mean(temp_surv)

        avg_vel = np.asarray(avg_vel)

        plot_ini = plotter_setup()
        fig, ax = plt.subplots()
        colour_axes = ax.scatter(in_pop, avg_vel, edgecolors='black', c = np.log10(avg_surv), zorder = 3)
        plt.colorbar(colour_axes, ax=ax, label = r'$\log_{10}\langle t_{ej}\rangle$ [Myr]')
        ax.set_ylabel(r'$\langle v_{ej} \rangle$ [km/s]')
        ax.axhline(676, color = 'black', linestyle = ':')
        ax.text(15, 690, r'$v_{esc, MW}$')
        ax.set_ylim(0, 1600)
        plot_ini.tickers_pop(ax, in_pop)
        plt.savefig('figures/ejection_stats/mean_vej.pdf', dpi=300, bbox_inches='tight')

        MW_code = MWpotentialBovy2015()
        vesc = (np.sqrt(2)*MW_code.circular_velocity(0.1 | units.parsec) + np.sqrt(2*constants.G*(4e6 | units.MSun)/(0.1 | units.parsec))).value_in(units.kms)
        
        fig, ax = plt.subplots()
        ax.set_title('Ejection Velocity Histogram')
        n, bins, patches = ax.hist(data.vesc, 20, histtype = 'step', color='black')
        n, bins, patches = ax.hist(data.vesc, 20, color='black', alpha = 0.3)
        
        ax.axvline(vesc, linestyle = ':', color = 'black')
        ax.text(vesc*(1+0.05), 0.9*max(n), r'$v_{esc, MW}$', horizontalalignment = 'center', rotation = 270)
        ax.set_xlabel(r'$v_{ejec}$ [km s$^{-1}$]')
        ax.set_ylabel(r'Occurence')
        ax.set_title('Ejection Velocity Histogram for All Simulations')
        plot_ini.tickers(ax, 'plot')
        plt.savefig('figures/ejection_stats/vejection_histogram.pdf', dpi = 300, bbox_inches='tight')

class event_tracker(object):
    """
    Class to take stats of ejection vs. mergers
    """
    def __init__(self):

        plot_ini = plotter_setup()
        chaos_data = glob.glob('data/Hermite/no_addition/chaotic_simulation/*')
        chaos_data = natsort.natsorted(chaos_data)

        init_pop = []
        merger = []
        for file_ in range(len(chaos_data)):
            with open(chaos_data[file_], 'rb') as input_file:
                data = pkl.load(input_file)
                init_pop.append(len(data.iloc[0][8])+1)
                merger.append(data.iloc[0][10])

        in_pop = np.unique(init_pop)
        frac_merge = np.empty(len(in_pop))
        iter = -1
        for pop_ in in_pop:
            iter +=1
            indices = np.where((init_pop == pop_))[0]
            temp_frac = [merger[i] for i in indices]
            frac_merge[iter] = np.mean(temp_frac)

        chaos_data_GRX = glob.glob('data/Hermite/no_addition/chaotic_simulation/*')
        chaos_data_GRX = natsort.natsorted(chaos_data)

        init_pop_GRX = []
        merger_GRX = []
        for file_ in range(len(chaos_data_GRX)):
            with open(chaos_data_GRX[file_], 'rb') as input_file:
                data = pkl.load(input_file)
                init_pop_GRX.append(len(data.iloc[0][8])+1)
                merger_GRX.append(data.iloc[0][10])

        in_pop_GRX = np.unique(init_pop_GRX)
        frac_merge_GRX = np.empty(len(in_pop_GRX))
        iter = -1
        for pop_ in in_pop_GRX:
            iter +=1
            indices = np.where((init_pop_GRX == pop_))[0]
            temp_frac = [merger_GRX[i] for i in indices]
            frac_merge_GRX[iter] = np.mean(temp_frac)

        fig, ax = plt.subplots()
        ax.set_ylabel(r'$N_{merge}/N_{sim}$')
        ax.set_ylim(0,1.05)
        ax.scatter(in_pop, frac_merge, color = 'red', edgecolors = 'black', label = 'Hermite')
        ax.scatter(in_pop_GRX, frac_merge_GRX, color = 'red', edgecolors = 'black', label = 'GRX')
        plot_ini.tickers_pop(ax, in_pop)
        plt.savefig('figures/ejection_stats/SMBH_merge_fraction.pdf', dpi=300, bbox_inches='tight')

#cst = KE_PE_plotters()
#cst.KEPE_plotter()

#cst = vejection()
#cst.vejec_histogram()
#cst.vejec_mean_plotter()
#cst = event_tracker()